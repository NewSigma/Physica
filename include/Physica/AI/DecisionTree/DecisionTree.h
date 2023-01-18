/*
 * Copyright 2023 WeiBo He.
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
#pragma once

#include <forward_list>
#include <set>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::AI {
    enum class DecisionTreeType {
        Regression,
        Classify
    };

    template<class ScalarType, DecisionTreeType Type> class RandomForest;
    /*
     * Reference:
     * [1] 周志华. 机器学习[M]. 清华大学出版社, 2016.73-
     */
    template<class ScalarType, DecisionTreeType Type>
    class DecisionTree {
        using MatrixType = Core::DenseMatrix<ScalarType>;
        using VectorType = Core::Vector<ScalarType>;
    public:
        struct Dataset {
            MatrixType features;
            VectorType labels;
            Utils::Array<bool> isFeatureContinuous;
        };

        using LossFunctor = ScalarType (*)(const Dataset&, const std::forward_list<size_t>&);
    private:
        size_t featureId;
        VectorType splitPoints;
        Utils::Array<DecisionTree> subTrees;
    public:
        DecisionTree(const Dataset& dataset);
        DecisionTree(const DecisionTree&) = default;
        DecisionTree(DecisionTree&&) noexcept = default;
        ~DecisionTree() = default;
        /* Operators */
        DecisionTree& operator=(DecisionTree tree) noexcept;
        /* Operations */
        [[nodiscard]] ScalarType predict(const VectorType& features) const;
        void swap(DecisionTree& tree) noexcept;
        /* Static members */
        static DecisionTree train(const Dataset& dataset);
    private:
        DecisionTree(size_t featureId_, VectorType splitPoints_, Utils::Array<DecisionTree> subTrees_);
        /* Getters */
        [[nodiscard]] bool isClassifyNode() const noexcept { return !isLeafNode() && splitPoints.getLength() != 1; }
        [[nodiscard]] bool isRegressionNode() const noexcept { return !isLeafNode() && splitPoints.getLength() == 1; }
        [[nodiscard]] bool isLeafNode() const noexcept { return subTrees.getLength() == 0; }
        /* Static members */
        static std::forward_list<size_t> makeInitialFeatures(size_t numFeature);
        static ScalarType checkStopRecursion(const Dataset& dataset,
                                             const std::forward_list<size_t>& availableSample,
                                             const std::forward_list<size_t>& availableFeature);
        template<class TrainFunctor>
        static DecisionTree doRecursion(const Dataset& dataset,
                                        const std::forward_list<size_t>& availableSample,
                                        const std::forward_list<size_t>& availableFeature,
                                        size_t featureId,
                                        VectorType splitPoints,
                                        TrainFunctor functor);
        static DecisionTree train(const Dataset& dataset,
                                  std::forward_list<size_t> availableSample,
                                  std::forward_list<size_t> availableFeature);
        static ScalarType findCommonLabel(const VectorType& labels, const std::forward_list<size_t>& availableSample);
        static ScalarType makeAverageLabel(const VectorType& labels, const std::forward_list<size_t>& availableSample);
        static ScalarType giniIndex(const Dataset& dataset, const std::forward_list<size_t>& availableSample);
        static ScalarType mse(const Dataset& dataset, const std::forward_list<size_t>& availableSample);
        static std::pair<size_t, VectorType> selectOptimalFeature(const Dataset& dataset,
                                                                  const std::forward_list<size_t>& availableSample,
                                                                  const std::forward_list<size_t>& availableFeature,
                                                                  LossFunctor functor);
        static inline LossFunctor getLossFunctor();
        /* Friends */
        friend class RandomForest<ScalarType, Type>;
    };

    template<class ScalarType, DecisionTreeType Type>
    inline typename DecisionTree<ScalarType, Type>::LossFunctor DecisionTree<ScalarType, Type>::getLossFunctor() {
        constexpr bool isClassifyTree = Type == DecisionTreeType::Classify;
        return isClassifyTree ? giniIndex : mse;
    }
}

#include "DecisionTreeImpl.h"
