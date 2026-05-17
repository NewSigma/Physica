//===- LifetimeMove.cpp - Narrowing lifetimes -----------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// The LifetimeMovePass identifies the precise lifetime range of allocas and
// repositions lifetime markers to stricter positions.
//
//===----------------------------------------------------------------------===//

#include "Physica/Transforms/LifetimeMovePass.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/Analysis/CaptureTracking.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/PostDominators.h"
#include "llvm/Analysis/PtrUseVisitor.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/Transforms/Coroutines/CoroInstr.h"

using namespace llvm;

#define DEBUG_TYPE "lifetime-move"

namespace {
    class LifetimeMover : public PtrUseVisitor<LifetimeMover> {
        using This = LifetimeMover;
        using Base = PtrUseVisitor<LifetimeMover>;

        const DominatorTree& DT;
        const LoopInfo& LI;

        SmallVector<AllocaInst*, 4> Allocas;
        // Critical points are instructions where the crossing of a variable's
        // lifetime makes a difference. We attempt to rise lifetime.end
        // before critical points and sink lifetime.start after them.
        SmallVector<Instruction*, 4> CriticalPoints;

        SmallVector<Instruction*, 2> LifetimeStarts;
        SmallVector<Instruction*, 2> LifetimeEnds;
        SmallPtrSet<BasicBlock*, 2> LifetimeStartBBs;

        SmallVector<Instruction*, 8> DirectUsers;
        SmallPtrSet<BasicBlock*, 2> DirectUserBBs;
    public:
        LifetimeMover(Function& F, const DominatorTree& DT, const LoopInfo& LI);
        LifetimeMover(const This&) = delete;

        bool run();

        void visitInstruction(Instruction& I);

        void visitAddrSpaceCastInst(AddrSpaceCastInst& ASC);
        void visitBitCastInst(BitCastInst& BC);
        void visitCallBase(CallBase& CB);
        void visitGetElementPtrInst(GetElementPtrInst& GEPI);
        void visitInsertElementInst(InsertElementInst& I);
        void visitInsertValueInst(InsertValueInst& I);
        void visitIntrinsicInst(IntrinsicInst& II);
        void visitPHINode(PHINode& I);
        void visitPtrToIntInst(PtrToIntInst& I);
        void visitMemIntrinsic(MemIntrinsic& I);
        void visitSelectInst(SelectInst& I);
        void visitStoreInst(StoreInst& SI);
    private:
        void addIndirectUser(Instruction& I);
        bool sinkLifetimeStartMarkers(AllocaInst* AI);
        bool riseLifetimeEndMarkers();
        void reset();
    };
} // namespace

LifetimeMover::LifetimeMover(Function& F, const DominatorTree& DT, const LoopInfo& LI)
        : Base(F.getDataLayout()), DT(DT), LI(LI) {
    for (Instruction& I : instructions(F)) {
        if (auto* AI = dyn_cast<AllocaInst>(&I)) {
            Allocas.push_back(AI);
            continue;
        }

        auto* II = dyn_cast<IntrinsicInst>(&I);
        if (II && II->isAssumeLikeIntrinsic())
            continue;

        if (isa<CallInst, InvokeInst, AnyCoroSuspendInst>(I))
            CriticalPoints.push_back(&I);
    }
}

bool LifetimeMover::run() {
    bool Changed = false;
    for (auto* AI : Allocas) {
        reset();
        Base::visitPtr(*AI);

        LLVM_DEBUG(AI->dump());
        if (!LifetimeStarts.empty()) {
            LLVM_DEBUG(dbgs() << "Try sink...\t");
            Changed |= sinkLifetimeStartMarkers(AI);
        }

        if (!LifetimeEnds.empty()) {
            LLVM_DEBUG(dbgs() << "Try rise...\t");
            // Do not move lifetime.end if alloca escapes
            if (PI.isEscaped())
                LLVM_DEBUG(dbgs() << "Escaped\n");
            else
                Changed |= riseLifetimeEndMarkers();
        }
    }
    return Changed;
}

void LifetimeMover::visitInstruction(Instruction& I) {
    DirectUsers.push_back(&I);
    DirectUserBBs.insert(I.getParent());
}

void LifetimeMover::visitAddrSpaceCastInst(AddrSpaceCastInst& ASC) {
    Base::visitAddrSpaceCastInst(ASC);
    addIndirectUser(ASC);
}

void LifetimeMover::visitBitCastInst(BitCastInst& BC) {
    Base::visitBitCastInst(BC);
    addIndirectUser(BC);
}

void LifetimeMover::visitCallBase(CallBase& CB) {
    switch (CB.getIntrinsicID()) {
    case Intrinsic::coro_await_suspend_bool:
    case Intrinsic::coro_await_suspend_handle:
    case Intrinsic::coro_await_suspend_void: {
        auto* AS = cast<CoroAwaitSuspendInst>(&CB);
        if (AS->getWrapperFunction()->arg_begin()->hasNoCaptureAttr()) {
            visitInstruction(CB);
            return;
        }
        break;
    }
    default:
    }

    for (unsigned Op = 0, OpCount = CB.arg_size(); Op < OpCount; ++Op)
        if (U->get() == CB.getArgOperand(Op) && !CB.doesNotCapture(Op))
            PI.setEscaped(&CB);
    InstVisitor<This>::visitCallBase(CB);
}

void LifetimeMover::visitGetElementPtrInst(GetElementPtrInst& GEPI) {
    Base::visitGetElementPtrInst(GEPI);
    addIndirectUser(GEPI);
}

void LifetimeMover::visitInsertElementInst(InsertElementInst& I) {
    enqueueUsers(I);
    addIndirectUser(I);
}

void LifetimeMover::visitInsertValueInst(InsertValueInst& I) {
    enqueueUsers(I);
    addIndirectUser(I);
}

void LifetimeMover::visitPHINode(PHINode& I) {
    enqueueUsers(I);
    addIndirectUser(I);
}

void LifetimeMover::visitIntrinsicInst(IntrinsicInst& II) {
    // lifetime markers are not actual uses
    switch (II.getIntrinsicID()) {
    case Intrinsic::lifetime_start:
        LifetimeStarts.push_back(&II);
        LifetimeStartBBs.insert(II.getParent());
        return;
    case Intrinsic::lifetime_end:
        LifetimeEnds.push_back(&II);
        return;
    default:
        Base::visitIntrinsicInst(II);
    }
}

void LifetimeMover::visitPtrToIntInst(PtrToIntInst& I) {
    Base::visitPtrToIntInst(I);
    visitInstruction(I);
}

void LifetimeMover::visitMemIntrinsic(MemIntrinsic& I) {
    visitInstruction(I);
}

void LifetimeMover::visitSelectInst(SelectInst& I) {
    enqueueUsers(I);
    addIndirectUser(I);
}

void LifetimeMover::visitStoreInst(StoreInst& SI) {
    if (SI.getPointerOperand() == U->get()) {
        InstVisitor<This>::visitStoreInst(SI);
        return;
    }

    PI.setEscaped(&SI);
    InstVisitor<This>::visitStoreInst(SI);
}
// TODO: Maybe we can move indirect users as well
void LifetimeMover::addIndirectUser(Instruction& I) {
    visitInstruction(I);
}
/// For each local variable that all of its user are dominated by one of the
/// critical point, we sink their lifetime.start markers to the place where
/// after the critical point. Doing so minimizes the lifetime of each variable.
bool LifetimeMover::sinkLifetimeStartMarkers(AllocaInst* AI) {
    auto Update = [this](Instruction* Old, Instruction* New) {
        if (Old == New)
            return Old;
        // Reject the new proposal if it lengthens lifetime
        if (DT.dominates(New, Old))
            return Old;

        if (LI.getLoopFor(New->getParent()))
            return Old;

        bool DomAll = llvm::all_of(DirectUserBBs, [this, New](BasicBlock* UserBB) {
            // Instruction level analysis if lifetime and users share a common BB
            BasicBlock* NewBB = New->getParent();
            if (UserBB == NewBB) {
                return llvm::all_of(DirectUsers, [New, UserBB](Instruction* I) {
                    return UserBB != I->getParent() || New->comesBefore(I);
                });
            }
            // Otherwise, BB level analysis is enough
            return DT.dominates(New, UserBB);
        });
        return DomAll ? New : Old;
    };

    // AllocaInst is a trivial critical point
    Instruction* DomPoint = AI;
    for (auto* P : CriticalPoints)
        DomPoint = Update(DomPoint, P);

    for (auto* P : DirectUsers)
        if (auto* N = P->getPrevNode())
            DomPoint = Update(DomPoint, N);

    // Sink lifetime.start markers to dominate block when they are
    // only used outside the region.
    if (DomPoint != AI) {
        // If existing position is better, do nothing
        for (auto* P : LifetimeStarts) {
            if (P == Update(DomPoint, P)) {
                LLVM_DEBUG(dbgs() << "Optimal\n");
                return false;
            }
        }

        // All the outsided lifetime.start markers are no longer necessary.
        auto* NewStart = LifetimeStarts[0]->clone();
        bool AnyErase = false;
        for (auto* I : LifetimeStarts) {
            if (DT.dominates(DomPoint, I) || LI.getLoopFor(I->getParent()))
                continue;

            bool Restart = llvm::any_of(LifetimeEnds, [this, I](Instruction* End) {
                return isPotentiallyReachable(End, I, &LifetimeStartBBs, &DT, &LI);
            });

            if (!Restart) {
                LifetimeStartBBs.erase(I->getParent());
                I->eraseFromParent();
                AnyErase = true;
            }
        }

        if (AnyErase) {
            LLVM_DEBUG(dbgs() << "Success: " << *DomPoint << '\n');
            if (DomPoint->isTerminator())
                NewStart->insertBefore(
                        cast<InvokeInst>(DomPoint)->getNormalDest()->getFirstNonPHIIt());
            else
                NewStart->insertAfter(DomPoint->getIterator());
        }
        else {
            // It is not beneficial if we failed to remove any lifetime markers
            LLVM_DEBUG(dbgs() << "No erase\n");
            NewStart->deleteValue();
        }
        return AnyErase;
    }
    LLVM_DEBUG(dbgs() << "Missing DomPoint\n");
    return false;
}
// Find the critical point that is dominated by all users of alloca,
// we will rise lifetime.end markers before the critical point.
bool LifetimeMover::riseLifetimeEndMarkers() {
    auto Update = [this](Instruction* Old, Instruction* New) {
        if (Old == New)
            return Old;

        if (Old != nullptr && DT.dominates(Old, New))
            return Old;

        if (LI.getLoopFor(New->getParent()))
            return Old;

        bool DomAll = llvm::all_of(DirectUserBBs, [this, New](BasicBlock* UserBB) {
            BasicBlock* NewBB = New->getParent();
            if (UserBB == NewBB) {
                return llvm::all_of(DirectUsers, [New, UserBB](Instruction* I) {
                    return UserBB != I->getParent() || I->comesBefore(New);
                });
            }

            if (auto* L = LI.getLoopFor(UserBB)) {
                SmallVector<BasicBlock*, 2> EBs;
                L->getOutermostLoop()->getExitingBlocks(EBs);
                return llvm::all_of(EBs, [this, NewBB](BasicBlock* EB) {
                    return DT.dominates(EB, NewBB);
                });
            }

            return DT.dominates(UserBB, NewBB);
        });
        return DomAll ? New : Old;
    };

    Instruction* DomPoint = nullptr;
    for (auto* P : CriticalPoints)
        DomPoint = Update(DomPoint, P);

    if (DomPoint != nullptr) {
        for (auto* P : LifetimeEnds) {
            if (P == Update(DomPoint, P)) {
                LLVM_DEBUG(dbgs() << "Optimal\n");
                return false;
            }
        }

        auto* NewEnd = LifetimeEnds[0]->clone();
        bool AnyErase = false;
        for (auto* I : LifetimeEnds) {
            if (!LI.getLoopFor(I->getParent())) {
                I->eraseFromParent();
                AnyErase = true;
            }
        }

        if (AnyErase) {
            LLVM_DEBUG(dbgs() << "Success: " << *DomPoint << '\n');
            NewEnd->insertBefore(DomPoint->getIterator());
        }
        else {
            LLVM_DEBUG(dbgs() << "No erase\n");
            NewEnd->deleteValue();
        }
        return AnyErase;
    }
    LLVM_DEBUG(dbgs() << "Missing DomPoint\n");
    return false;
}

void LifetimeMover::reset() {
    PI.reset();
    Worklist.clear();
    VisitedUses.clear();

    LifetimeStarts.clear();
    LifetimeEnds.clear();
    LifetimeStartBBs.clear();

    DirectUsers.clear();
    DirectUserBBs.clear();
}

PreservedAnalyses LifetimeMovePass::run(Function& F, FunctionAnalysisManager& AM) {
    const DominatorTree& DT = AM.getResult<DominatorTreeAnalysis>(F);
    const LoopInfo& LI = AM.getResult<LoopAnalysis>(F);
    LifetimeMover Mover(F, DT, LI);
    if (!Mover.run())
        return PreservedAnalyses::all();

    PreservedAnalyses PA;
    PA.preserveSet<CFGAnalyses>();
    return PA;
}
