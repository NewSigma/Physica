/*
 * Copyright 2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "Physica/Transforms/LifetimeMovePass.h"

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
    __attribute__((visibility("default"))) llvm::PassPluginLibraryInfo llvmGetPassPluginInfo() {
        return llvm::PassPluginLibraryInfo{
            .APIVersion = LLVM_PLUGIN_API_VERSION,
            .PluginName = "PhysicaTransforms",
            .PluginVersion = LLVM_VERSION_STRING,
            .RegisterPassBuilderCallbacks = registerPassBuilderCallbacks,
        };
    }
}
