/*
 * Copyright 2023-2024 Weibo He.
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

namespace Physica::Core {
    enum class DecisionTreeType {
        Regression,
        Classify
    };

    template<Scalar T, DecisionTreeType Type> class RandomForest;
    /*
     * Reference:
     * [1] 周志华. 机器学习[M]. 清华大学出版社, 2016.73-
     */
    template<Scalar T, DecisionTreeType Type>
    class DecisionTree {
        using MatrixType = DenseMatrix<T>;
        using VectorType = VectorND<T>;
    public:
        struct Dataset {
            MatrixType features;
            VectorType labels;
            Array<bool> isFeatureContinuous;
        };

        using LossFunctor = T (*)(const Dataset&, const Array<size_t>&);
    private:
        enum NodeType {
            Regression,
            Classify
        };

        size_t featureId;
        T splitPoint;
        Array<DecisionTree> subTrees;
        NodeType nodeType;
    public:
        DecisionTree(const Dataset& dataset);
        DecisionTree(const DecisionTree&) = default;
        DecisionTree(DecisionTree&&) noexcept = default;
        ~DecisionTree() = default;
        /* Operators */
        DecisionTree& operator=(DecisionTree tree) noexcept;
        /* Operations */
        [[nodiscard]] T predict(const VectorType& features) const;
        void swap(DecisionTree& __restrict tree) noexcept;
        /* Static members */
        static DecisionTree train(const Dataset& dataset);
    private:
        DecisionTree(size_t featureId_, T splitPoint_, Array<DecisionTree> subTrees_, NodeType nodeType_);
        /* Getters */
        [[nodiscard]] bool isClassifyNode() const noexcept { return !isLeafNode() && nodeType == NodeType::Classify; }
        [[nodiscard]] bool isRegressionNode() const noexcept { return !isLeafNode() && nodeType == NodeType::Regression; }
        [[nodiscard]] bool isLeafNode() const noexcept { return subTrees.getLength() == 0; }
        /* Static members */
        static std::forward_list<size_t> makeInitialFeatures(size_t numFeature);
        static T checkStopRecursion(const Dataset& dataset,
                                             const std::forward_list<size_t>& availableSample,
                                             const std::forward_list<size_t>& availableFeature);
        template<class TrainFunctor>
        static DecisionTree doRecursion(const Dataset& dataset,
                                        const std::forward_list<size_t>& availableSample,
                                        const std::forward_list<size_t>& availableFeature,
                                        size_t featureId,
                                        T splitPoint,
                                        TrainFunctor functor);
        static DecisionTree train(const Dataset& dataset,
                                  std::forward_list<size_t> availableSample,
                                  std::forward_list<size_t> availableFeature);
        static T findCommonLabel(const VectorType& labels, const std::forward_list<size_t>& availableSample);
        static T makeAverageLabel(const VectorType& labels, const std::forward_list<size_t>& availableSample);
        static T giniIndex(const Dataset& dataset, const Array<size_t>& availableSample);
        static T mse(const Dataset& dataset, const Array<size_t>& availableSample);
        static std::pair<size_t, T> selectOptimalFeature(const Dataset& dataset,
                                                                  const std::forward_list<size_t>& availableSample,
                                                                  const std::forward_list<size_t>& availableFeature,
                                                                  LossFunctor functor);
        static inline LossFunctor getLossFunctor();
        /* Friends */
        friend class RandomForest<T, Type>;
    };

    template<Scalar T, DecisionTreeType Type>
    inline typename DecisionTree<T, Type>::LossFunctor DecisionTree<T, Type>::getLossFunctor() {
        constexpr bool isClassifyTree = Type == DecisionTreeType::Classify;
        return isClassifyTree ? giniIndex : mse;
    }
}

#include "DecisionTreeImpl.h"
