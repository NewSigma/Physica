/*
 * Copyright 2023 Weibo He.
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

namespace Physica::Core {
    /*
     * Reference:
     * [1] Leo Breiman, Random forests[J]. Machine Learning, 45, 5–32, 2001
     */

    template<class ScalarType, DecisionTreeType Type>
    class RandomForest {
        using TreeType = DecisionTree<ScalarType, Type>;
        using VectorType = typename TreeType::VectorType;
        using MatrixType = typename TreeType::MatrixType;
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
        void swap(RandomForest& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumTree() const noexcept { return trees.getLength(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] static std::pair<RandomForest, ScalarType> train(unsigned int numTree, const Dataset& dataset, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static ScalarType makeFeatureImportance(
                size_t featureId,
                unsigned int numForest,
                unsigned int numTree,
                Dataset dataset,
                RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] std::forward_list<size_t> selectImportantFeature(
                unsigned int numForest,
                unsigned int numTree,
                Dataset dataset,
                RandomGenerator& gen);
    private:
        RandomForest(Utils::Array<TreeType> trees_);
        template<class RandomGenerator>
        static TreeType trainTree(const Dataset& dataset,
                                  std::forward_list<size_t> availableSample,
                                  std::forward_list<size_t> availableFeature,
                                  RandomGenerator& gen);
        template<class RandomGenerator>
        static std::pair<std::forward_list<size_t>, std::forward_list<size_t>> randTrainTestSet(size_t numSample, RandomGenerator& gen);
        static void testTree(ScalarType prediction, ScalarType label, size_t sampleId, VectorType& predictions, Utils::Array<size_t>& numTestSample);
        static ScalarType makeTestError(const VectorType& predictions, const VectorType& labels, const Utils::Array<size_t>& numTestSample);
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
                toNextMean(mean, i, trees[i].predict(features));
            return mean;
        }
    }

    template<class ScalarType, DecisionTreeType Type>
    void RandomForest<ScalarType, Type>::swap(RandomForest& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        trees.swap(obj.trees);
    }

    template<class ScalarType, DecisionTreeType Type>
    template<class RandomGenerator>
    std::pair<RandomForest<ScalarType, Type>, ScalarType> RandomForest<ScalarType, Type>::train(
            unsigned int numTree,
            const Dataset& dataset,
            RandomGenerator& gen) {
        const size_t numSample = dataset.features.getRow();
        const auto availableFeature = TreeType::makeInitialFeatures(dataset.features.getColumn());
        Utils::Array<TreeType> trees{};
        trees.reserve(numTree);
        VectorType predictions(numSample, 0);
        Utils::Array<size_t> numTestSample(numSample, 0);
        for (size_t i = 0; i < numTree; ++i) {
            auto pair = randTrainTestSet(numSample, gen);
            auto& trainSample = pair.first;
            auto& testSample = pair.second;

            auto tree = trainTree(dataset, std::move(trainSample), availableFeature, gen);
            for (size_t sample : testSample)
                testTree(tree.predict(dataset.features.row(sample)), dataset.labels[sample], sample, predictions, numTestSample);
            trees.append(std::move(tree));
        }
        return {RandomForest(std::move(trees)), makeTestError(predictions, dataset.labels, numTestSample)};
    }
    /*
     * Reference
     * [1] Kursa M B, Rudnicki W R. Feature selection with the Boruta package[J]. J Stat Softw, 2010, 36(11): 1-13.
     * [2] Manual On Setting Up, Using, And Understanding Random Forests V3.1 (https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf)
     */
    template<class ScalarType, DecisionTreeType Type>
    template<class RandomGenerator>
    ScalarType RandomForest<ScalarType, Type>::makeFeatureImportance(
            size_t featureId,
            unsigned int numForest,
            unsigned int numTree,
            Dataset dataset,
            RandomGenerator& gen) {
        const size_t numSample = dataset.features.getRow();
        const size_t numFeature = dataset.features.getColumn();
        assert(featureId < numFeature * 2);
        const auto availableFeature = TreeType::makeInitialFeatures(numFeature);
        auto& features = dataset.features;
        std::uniform_int_distribution<size_t> dist(0, numSample - 1);
        /* Extend features */ {
            features.resize(numSample, numFeature * 2);
            dataset.isFeatureContinuous.resize(features.getColumn());
            for (size_t i = 0; i < numFeature; ++i)
                dataset.isFeatureContinuous[i + numFeature] = dataset.isFeatureContinuous[i];
        }
        ScalarType mean = 0;
        for (size_t i = 0; i < numForest; ++i) {
            for (size_t i = 0; i < numFeature; ++i)
                for (size_t j = 0; j < numSample; ++j)
                    features(j, i + numFeature) = features(dist(gen), i);

            VectorType predictions(numSample, 0);
            Utils::Array<size_t> numTestSample(numSample, 0);
            for (size_t _ = 0; _ < numTree; ++_) {
                auto pair = randTrainTestSet(numSample, gen);
                auto& trainSample = pair.first;
                auto& testSample = pair.second;

                auto tree = trainTree(dataset, std::move(trainSample), availableFeature, gen);
                for (size_t sample : testSample) {
                    VectorType feature = features.row(sample);
                    feature[featureId] = features(dist(gen), featureId);
                    testTree(tree.predict(feature), dataset.labels[sample], sample, predictions, numTestSample);
                }
            }
            const ScalarType newError = makeTestError(predictions, dataset.labels, numTestSample);
            const ScalarType deltaError = newError - train(numTree, dataset, gen).second;
            Core::toNextMean(mean, i, deltaError);
        }
        return mean;
    }

    template<class ScalarType, DecisionTreeType Type>
    template<class RandomGenerator>
    std::forward_list<size_t> RandomForest<ScalarType, Type>::selectImportantFeature(
            unsigned int numForest,
            unsigned int numTree,
            Dataset dataset,
            RandomGenerator& gen) {
        std::forward_list<size_t> result{};
        for (size_t i = 0; i < dataset.isFeatureContinuous.getLength(); ++i)
            result.push_front(i);
        
        const size_t numSample = dataset.labels.getLength();
        while (true) {
            const size_t numFeature = std::distance(result.begin(), result.end());
            MatrixType features(numSample, numFeature);
            Utils::Array<bool> isFeatureContinuous(numFeature);
            {
                size_t index = 0;
                for (auto ite = result.begin(); ite != result.end(); ++ite) {
                    features.col(index) = dataset.features.col(*ite).asVector();
                    isFeatureContinuous[index] = dataset.isFeatureContinuous[*ite];
                }
            }
            auto importance = VectorType(numFeature * 2);
            for (size_t i = 0; i < importance.getLength(); ++i)
                importance[i] = makeFeatureImportance(i, numForest, numTree, {features, dataset.labels, isFeatureContinuous}, gen);
            
            const auto shadow_importance = importance.tail(numFeature);
            const ScalarType upper_bound = mean(shadow_importance) + deviation(shadow_importance);
            auto ite = result.before_begin();
            for (size_t i = numFeature - 1; i < numFeature; --i, ++ite) {
                if (importance[i] < upper_bound)
                    result.erase_after(ite);
            }
        }
        return result;
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
            return TreeType(dataset.features.getColumn(), criteria, {}, TreeType::NodeType::Classify);

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

    template<class ScalarType, DecisionTreeType Type>
    template<class RandomGenerator>
    std::pair<std::forward_list<size_t>, std::forward_list<size_t>>
    RandomForest<ScalarType, Type>::randTrainTestSet(size_t numSample, RandomGenerator& gen) {
        std::uniform_int_distribution<size_t> dist(0, numSample - 1);
        std::set<size_t> set{};
        for (size_t _ = 0; _ < numSample; ++_)
            set.insert(dist(gen));

        std::forward_list<size_t> trainSample{}, testSample{};
        for (size_t j = 0; j < numSample; ++j) {
            if (set.find(j) == set.end())
                testSample.push_front(j);
            else
                trainSample.push_front(j);
        }
        return std::make_pair(std::move(trainSample), std::move(testSample));
    }

    template<class ScalarType, DecisionTreeType Type>
    void RandomForest<ScalarType, Type>::testTree(
            ScalarType prediction,
            [[maybe_unused]] ScalarType label,
            size_t sampleId,
            VectorType& predictions,
            Utils::Array<size_t>& numTestSample) {
        if constexpr (Type == DecisionTreeType::Classify) {
            const bool isCorrect = prediction == label;
            predictions[sampleId] += ScalarType(isCorrect ? 1.0 : -1.0);
            numTestSample[sampleId] += 1;
        }
        else {
            toNextMean(predictions[sampleId], numTestSample[sampleId], prediction);
            numTestSample[sampleId] += 1;
        }
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType RandomForest<ScalarType, Type>::makeTestError(
            const VectorType& predictions,
            const VectorType& labels,
            [[maybe_unused]] const Utils::Array<size_t>& numTestSample) {
        const size_t numSample = predictions.getLength();
        if constexpr (Type == DecisionTreeType::Classify) {
            size_t count = 0;
            for (size_t i = 0; i < numSample; ++i) {
                const bool isPredictionCorrect = predictions[i].isPositive();
                count += isPredictionCorrect;
            }
            const ScalarType accuracy = ScalarType(count) / ScalarType(numSample);
            return accuracy;
        }
        else {
            ScalarType mse = 0;
            for (size_t i = 0; i < numSample; ++i)
                mse += square(predictions[i] - labels[i]);
            mse /= ScalarType(numSample);
            return mse;
        }
    }
}
