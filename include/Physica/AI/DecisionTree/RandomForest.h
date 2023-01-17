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

#include "DecisionTree.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::AI {
    /*
     * Reference:
     * [1] Leo Breiman, Random forests[J]. Machine Learning, 45, 5–32, 2001
     */

    template<class ScalarType, DecisionTreeType Type>
    class RandomForest {
        using TreeType = DecisionTree<ScalarType, Type>;
        using VectorType = typename TreeType::VectorType;
    public:
        using Dataset = typename TreeType::Dataset;
    private:
        Utils::Array<TreeType> trees;
    public:
        RandomForest(const RandomForest&) = default;
        RandomForest(RandomForest&&) noexcept = default;
        ~RandomForest() = default;
        /* Operators */
        RandomForest& operator=(RandomForest obj) noexcept;
        /* Operations */
        [[nodiscard]] ScalarType predict(const VectorType& features) const;
        void swap(RandomForest& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumTrees() const noexcept { return trees.getLength(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] static std::pair<RandomForest, ScalarType> train(unsigned int numTree, const Dataset& dataset, RandomGenerator& gen);
    private:
        RandomForest(Utils::Array<TreeType> trees_);
        template<class RandomGenerator>
        static TreeType trainTree(const Dataset& dataset,
                                  std::forward_list<size_t> availableSample,
                                  std::forward_list<size_t> availableFeature,
                                  RandomGenerator& gen);
    };

    template<class ScalarType, DecisionTreeType Type>
    RandomForest<ScalarType, Type>::RandomForest(Utils::Array<TreeType> trees_) : trees(std::move(trees_)) {}

    template<class ScalarType, DecisionTreeType Type>
    RandomForest<ScalarType, Type>& RandomForest<ScalarType, Type>::operator=(RandomForest obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType RandomForest<ScalarType, Type>::predict(const VectorType& features) const {
        if constexpr (Type == DecisionTreeType::Classify) {
            std::multiset<ScalarType> labels{};
            for (const auto& tree : trees)
                labels.insert(tree.predict(features));
            
            size_t maxCount = 0;
            ScalarType result = 0;
            for (auto ite = labels.begin(); ite != labels.end(); ++ite) {
                const auto count = labels.count(*ite);
                if (count > maxCount) {
                    maxCount = count;
                    result = *ite;
                }
            }
            return result;
        }
        else {
            ScalarType mean = 0;
            for (size_t i = 0; i < trees.getLength(); ++i)
                toNextMean(mean, trees[i].predict(features), i);
            return mean;
        }
    }

    template<class ScalarType, DecisionTreeType Type>
    void RandomForest<ScalarType, Type>::swap(RandomForest& obj) noexcept {
        trees.swap(obj.trees);
    }

    template<class ScalarType, DecisionTreeType Type>
    template<class RandomGenerator>
    std::pair<RandomForest<ScalarType, Type>, ScalarType> RandomForest<ScalarType, Type>::train(
            unsigned int numTree,
            const Dataset& dataset,
            RandomGenerator& gen) {
        std::forward_list<size_t> availableFeature{};
        const size_t numFeature = dataset.features.getColumn();
        for (size_t i = 0; i < numFeature; ++i)
            availableFeature.push_front(i);

        const size_t numSample = dataset.features.getRow();
        Utils::Array<TreeType> trees{};
        trees.reserve(numTree);
        ScalarType forestError = 0;
        std::uniform_int_distribution<size_t> dist(0, numSample - 1);
        for (size_t i = 0; i < numTree; ++i) {
            std::set<size_t> set{};
            for (size_t _ = 0; _ < numSample; ++_)
                set.insert(dist(gen));

            std::forward_list<size_t> trainSample{}, testSample{};
            size_t numTestSample = 0;
            for (size_t j = 0; j < numSample; ++j) {
                if (set.find(j) == set.end()) {
                    testSample.push_front(j);
                    numTestSample += 1;
                }
                else
                    trainSample.push_front(j);
            }

            auto tree = trainTree(dataset, std::move(trainSample), availableFeature, gen);
            ScalarType treeError;
            if constexpr (Type == DecisionTreeType::Classify) {
                size_t count = 0;
                for (size_t sample : testSample)
                    count += tree.predict(dataset.features.row(sample)) == dataset.labels[sample];
                treeError = ScalarType(count) / ScalarType(numTestSample);
            }
            else {
                treeError = 0;
                for (size_t sample : testSample)
                    treeError += square(tree.predict(dataset.features.row(sample)) - dataset.labels[sample]);
                treeError /= ScalarType(numTestSample);
            }
            Core::toNextMean(forestError, i, treeError);
            trees.append(std::move(tree));
        }
        return {RandomForest(std::move(trees)), forestError};
    }

    template<class ScalarType, DecisionTreeType Type>
    template<class RandomGenerator>
    typename RandomForest<ScalarType, Type>::TreeType RandomForest<ScalarType, Type>::trainTree(
            const Dataset& dataset,
            std::forward_list<size_t> availableSample,
            std::forward_list<size_t> availableFeature,
            RandomGenerator& gen) {
        std::forward_list<size_t> randomFeature = availableFeature;
        /* make randomFeature */ {
            const size_t numAvailableFeature = std::distance(availableFeature.cbegin(), availableFeature.cend());
            const size_t numRandFeature = double(ln(ScalarType(numAvailableFeature)) / M_LN2 + 1.0);
            const bool isResultValid = numAvailableFeature > 0;
            const bool shouldRemoveElem = numRandFeature != numAvailableFeature;
            if (isResultValid && shouldRemoveElem) {
                for (size_t i = 0; i < (numAvailableFeature - numRandFeature); ++i) {
                    std::uniform_int_distribution<size_t> dist(0, numAvailableFeature - i - 1);
                    auto ite = randomFeature.before_begin();
                    std::advance(ite, dist(gen));
                    randomFeature.erase_after(ite);
                }
            }
        }
        const ScalarType criteria = TreeType::checkStopRecursion(dataset, availableSample, randomFeature);
        const bool shouldStopRecursion = !std::isnan(double(criteria));
        if (shouldStopRecursion)
            return TreeType(dataset.features.getColumn(), {criteria}, {});

        auto pair = TreeType::selectOptimalFeature(dataset,
                                                   availableSample,
                                                   randomFeature,
                                                   TreeType::getLossFunctor());
        return TreeType::doRecursion(dataset,
                                     availableSample,
                                     availableFeature,
                                     pair.first,
                                     std::move(pair.second),
                                     [&gen](const Dataset& dataset,
                                            std::forward_list<size_t> availableSample,
                                            std::forward_list<size_t> availableFeature) {
                                        return trainTree(dataset, availableSample, availableFeature, gen);
                                     });
    }
}
