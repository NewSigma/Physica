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

namespace Physica::AI {
    template<class ScalarType, DecisionTreeType Type>
    DecisionTree<ScalarType, Type>::DecisionTree(const Dataset& dataset) : DecisionTree(train(dataset)) {}

    template<class ScalarType, DecisionTreeType Type>
    DecisionTree<ScalarType, Type>::DecisionTree(size_t featureId_, VectorType splitPoints_, Utils::Array<DecisionTree> subTrees_)
            : featureId(featureId_), splitPoints(std::move(splitPoints_)), subTrees(std::move(subTrees_)) {}

    template<class ScalarType, DecisionTreeType Type>
    DecisionTree<ScalarType, Type>& DecisionTree<ScalarType, Type>::operator=(DecisionTree tree) noexcept {
        swap(tree);
        return *this;
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType DecisionTree<ScalarType, Type>::predict(const VectorType& features) const {
        if (isLeafNode())
            return splitPoints[0];

        const ScalarType feature = features[featureId];
        if (isClassifyNode()) {
            for (size_t i = 0; i < splitPoints.getLength(); ++i) {
                if (feature == splitPoints[i])
                    return subTrees[i].predict(features);
            }
        }
        else {
            const size_t index = (feature < splitPoints[0]) ? 0U : 1U;
            return subTrees[index].predict(features);
        }
        return std::nan("");
    }

    template<class ScalarType, DecisionTreeType Type>
    void DecisionTree<ScalarType, Type>::swap(DecisionTree& tree) noexcept {
        std::swap(featureId, tree.featureId);
        splitPoints.swap(tree.splitPoints);
        subTrees.swap(tree.subTrees);
    }

    template<class ScalarType, DecisionTreeType Type>
    DecisionTree<ScalarType, Type> DecisionTree<ScalarType, Type>::train(const Dataset& dataset) {
        std::forward_list<size_t> availableSample;
        std::forward_list<size_t> availableFeature;
        const size_t numSample = dataset.features.getRow();
        const size_t numFeature = dataset.features.getColumn();
        for (size_t i = 0; i < numSample; ++i)
            availableSample.push_front(i);
        for (size_t i = 0; i < numFeature; ++i)
            availableFeature.push_front(i);
        return train(dataset, std::move(availableSample), std::move(availableFeature));
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType DecisionTree<ScalarType, Type>::checkStopRecursion(
            const Dataset& dataset,
            const std::forward_list<size_t>& availableSample,
            const std::forward_list<size_t>& availableFeature) {
        assert(!availableSample.empty());
        constexpr bool isClassifyTree = Type == DecisionTreeType::Classify;
        const auto& features = dataset.features;
        const auto& labels = dataset.labels;

        if constexpr (isClassifyTree) {
            const ScalarType label = labels[availableSample.front()];
            bool isAllLabelSame = true;
            for (size_t index : availableSample)
                isAllLabelSame &= label == labels[index];
            
            if (isAllLabelSame)
                return label;
        }

        const bool isFeatureEmpty = availableFeature.empty();
        if (isFeatureEmpty) {
            ScalarType prediction;
            if constexpr (isClassifyTree)
                prediction = findCommonLabel(labels, availableSample);
            else
                prediction = makeAverageLabel(labels, availableSample);
            return prediction;
        }

        for (auto featureId : availableFeature) {
            const ScalarType feature = features(availableSample.front(), featureId);
            for (size_t sample : availableSample) {
                const bool isSeparable = features(sample, featureId) != feature;
                if (isSeparable)
                    return std::nan("");
            }
        }

        if constexpr (Type == DecisionTreeType::Classify)
            return findCommonLabel(labels, availableSample);
        else
            return makeAverageLabel(labels, availableSample);
    }

    template<class ScalarType, DecisionTreeType Type>
    template<class TrainFunctor>
    DecisionTree<ScalarType, Type> DecisionTree<ScalarType, Type>::doRecursion(const Dataset& dataset,
                                                                               const std::forward_list<size_t>& availableSample,
                                                                               const std::forward_list<size_t>& availableFeature,
                                                                               size_t featureId,
                                                                               VectorType splitPoints,
                                                                               TrainFunctor functor) {
        std::forward_list<size_t> newAvailableFeature{};
        for (auto feature : availableFeature) {
            if (feature != featureId)
                newAvailableFeature.push_front(feature);
        }

        Utils::Array<DecisionTree> subTrees{};
        if (dataset.isFeatureContinuous[featureId]) {
            std::forward_list<size_t> list1{}, list2{};
            for (auto sample : availableSample) {
                if (dataset.features(sample, featureId) < splitPoints[0])
                    list1.push_front(sample);
                else
                    list2.push_front(sample);
            }

            subTrees.reserve(2U);
            subTrees.append(functor(dataset, std::move(list1), newAvailableFeature));
            subTrees.append(functor(dataset, std::move(list2), newAvailableFeature));
        }
        else {
            subTrees.reserve(splitPoints.getLength());
            for (auto point : splitPoints) {
                std::forward_list<size_t> newAvailableSample{};
                for (auto sample : availableSample) {
                    if (dataset.features(sample, featureId) == point)
                        newAvailableSample.push_front(sample);
                }
                auto tree = functor(dataset, std::move(newAvailableSample), std::move(newAvailableFeature));
                subTrees.append(std::move(tree));
            }
        }
        return DecisionTree(featureId, std::move(splitPoints), std::move(subTrees));
    }

    template<class ScalarType, DecisionTreeType Type>
    DecisionTree<ScalarType, Type> DecisionTree<ScalarType, Type>::train(
            const Dataset& dataset,
            std::forward_list<size_t> availableSample,
            std::forward_list<size_t> availableFeature) {
        using TrainFunctor = DecisionTree (*)(const Dataset& dataset, std::forward_list<size_t>, std::forward_list<size_t>);
        const size_t numFeature = dataset.features.getColumn();
        const ScalarType criteria = checkStopRecursion(dataset, availableSample, availableFeature);
        const bool shouldStopRecursion = !std::isnan(double(criteria));
        if (shouldStopRecursion)
            return DecisionTree(numFeature, {criteria}, {});

        constexpr bool isClassifyTree = Type == DecisionTreeType::Classify;
        auto pair = selectOptimalFeature(dataset, availableSample, availableFeature, isClassifyTree ? giniIndex : mse);
        return doRecursion<TrainFunctor>(dataset, availableSample, availableFeature, pair.first, std::move(pair.second), train);
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType DecisionTree<ScalarType, Type>::giniIndex(const Dataset& dataset, const std::forward_list<size_t>& availableSample) {
        std::multiset<ScalarType> set{};
        for (size_t sample : availableSample)
            set.insert(dataset.labels[sample]);
        
        ScalarType result = 1;
        const ScalarType numSample = std::distance(availableSample.begin(), availableSample.end());
        for (auto ite = set.begin(); ite != set.end(); ++ite) {
            const ScalarType temp = set.count(*ite);
            result -= square(temp / numSample);
        }
        return result;
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType DecisionTree<ScalarType, Type>::mse(const Dataset& dataset, const std::forward_list<size_t>& availableSample) {
        ScalarType sum = 0;
        ScalarType sumSquare = 0;
        size_t count = 0;
        for (auto sample : availableSample) {
            const ScalarType label = dataset.labels[sample];
            sum += label;
            sumSquare += square(label);
            count += 1;
        }
        return ScalarType(count) * sumSquare - square(sum);
    }

    template<class ScalarType, DecisionTreeType Type>
    std::pair<size_t, typename DecisionTree<ScalarType, Type>::VectorType>
    DecisionTree<ScalarType, Type>::selectOptimalFeature(
            const Dataset& dataset,
            const std::forward_list<size_t>& availableSample,
            const std::forward_list<size_t>& availableFeature,
            LossFunctor functor) {
        assert(!availableFeature.empty());
        size_t optimalFeatureId = dataset.features.getColumn();
        VectorType splitPoints{};
        ScalarType minLoss = std::numeric_limits<ScalarType>::max();
        for (auto featureId : availableFeature) {
            VectorType featureVector{};
            {
                std::set<ScalarType> featureSet{};
                for (auto sample : availableSample)
                    featureSet.insert(dataset.features(sample, featureId));

                featureVector.resize(featureSet.size());
                size_t i = 0;
                for (auto ite = featureSet.begin(); ite != featureSet.end(); ++ite) {
                    featureVector[i] = *ite;
                    i += 1;
                }
            }

            ScalarType loss = 0;
            if (dataset.isFeatureContinuous[featureId]) {
                std::forward_list<size_t> list1{}, list2{};
                loss = std::numeric_limits<ScalarType>::max();
                ScalarType splitPoint = 0;
                for (size_t i = 0; i < featureVector.getLength() - 1; ++i) {
                    const ScalarType midpoint = (featureVector[i] + featureVector[i + 1]) * 0.5;
                    size_t weight1 = 0;
                    for (auto sample : availableSample) {
                        if (dataset.features(sample, featureId) < midpoint) {
                            list1.push_front(sample);
                            weight1 += 1;
                        }
                        else
                            list2.push_front(sample);
                    }

                    const size_t weight2 = featureVector.getLength() - weight1;
                    const ScalarType temp = functor(dataset, list1) * ScalarType(weight1) + functor(dataset, list2) * ScalarType(weight2);
                    if (temp < loss) {
                        loss = temp;
                        splitPoint = midpoint;
                    }
                    list1.clear();
                    list2.clear();
                }

                if (loss <= minLoss) {
                    minLoss = loss;
                    optimalFeatureId = featureId;
                    splitPoints = {splitPoint};
                }
            }
            else {
                std::forward_list<size_t> list{};
                loss = 0;
                for (auto value : featureVector) {
                    size_t weight = 0;
                    for (auto sample : availableSample) {
                        if (dataset.features(sample, featureId) == value) {
                            list.push_front(sample);
                            weight += 1;
                        }
                    }
                    loss += functor(dataset, list) * ScalarType(weight);
                    list.clear();
                }

                if (loss <= minLoss) {
                    minLoss = loss;
                    optimalFeatureId = featureId;
                    splitPoints.swap(featureVector);
                }
            }
        }
        return {optimalFeatureId, std::move(splitPoints)};
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType DecisionTree<ScalarType, Type>::findCommonLabel(
            const VectorType& labels,
            const std::forward_list<size_t>& availableSample) {
        assert(!availableSample.empty());
        std::multiset<ScalarType> set{};
        for (size_t sample : availableSample)
            set.insert(labels[sample]);
        
        size_t numCount = 0;
        ScalarType result = 0;
        for (auto ite = set.begin(); ite != set.end(); ++ite) {
            size_t temp = set.count(*ite);
            if (temp > numCount) {
                numCount = temp;
                result = *ite;
            }
        }
        return result;
    }

    template<class ScalarType, DecisionTreeType Type>
    ScalarType DecisionTree<ScalarType, Type>::makeAverageLabel(const VectorType& labels, const std::forward_list<size_t>& availableSample) {
        ScalarType result = 0;
        size_t count = 0;
        for (auto sample : availableSample) {
            result += labels[sample];
            count += 1;
        }
        result /= ScalarType(count);
        return result;
    }
}
