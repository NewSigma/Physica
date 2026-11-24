/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "LifetimeMovePass.h"

using namespace llvm;

namespace {
    void registerPassBuilderCallbacks(PassBuilder& PB) {
        PB.registerPeepholeEPCallback([](llvm::FunctionPassManager& PM, OptimizationLevel Level) {
            if (Level == OptimizationLevel::O2 || Level == OptimizationLevel::O3)
                PM.addPass(LifetimeMovePass());
        });
    }
}

extern "C" {
    [[gnu::visibility("default")]] llvm::PassPluginLibraryInfo llvmGetPassPluginInfo() {
        return llvm::PassPluginLibraryInfo{
            .APIVersion = LLVM_PLUGIN_API_VERSION,
            .PluginName = "PhysicaTransforms",
            .PluginVersion = LLVM_VERSION_STRING,
            .RegisterPassBuilderCallbacks = registerPassBuilderCallbacks,
        };
    }
}
