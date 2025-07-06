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

namespace Physica {
    template<Scalar T, DecisionTreeType Type>
    DecisionTree<T, Type>::DecisionTree(const Dataset& dataset) : DecisionTree(train(dataset)) {}

    template<Scalar T, DecisionTreeType Type>
    DecisionTree<T, Type>::DecisionTree(size_t featureId_, T splitPoint_, Array<DecisionTree> subTrees_, NodeType nodeType_)
            : featureId(featureId_)
            , splitPoint(std::move(splitPoint_))
            , subTrees(std::move(subTrees_))
            , nodeType(nodeType_) {}

    template<Scalar T, DecisionTreeType Type>
    DecisionTree<T, Type>& DecisionTree<T, Type>::operator=(DecisionTree tree) noexcept {
        swap(tree);
        return *this;
    }

    template<Scalar T, DecisionTreeType Type>
    T DecisionTree<T, Type>::predict(const VectorType& features) const {
        if (isLeafNode())
            return splitPoint;

        const T feature = features[featureId];
        size_t index;
        if (isClassifyNode())
            index = (feature == splitPoint) ? 0U : 1U;
        else
            index = (feature < splitPoint) ? 0U : 1U;
        return subTrees[index].predict(features);
    }

    template<Scalar T, DecisionTreeType Type>
    void DecisionTree<T, Type>::swap(DecisionTree& __restrict tree) noexcept {
        assert(this != &tree && "[Error]: Self swap is likely a bug");
        std::swap(featureId, tree.featureId);
        splitPoint.swap(tree.splitPoint);
        subTrees.swap(tree.subTrees);
        std::swap(nodeType, tree.nodeType);
    }

    template<Scalar T, DecisionTreeType Type>
    DecisionTree<T, Type> DecisionTree<T, Type>::train(const Dataset& dataset) {
        std::forward_list<size_t> availableSample;
        const size_t numSample = dataset.features.getRow();
        for (size_t i = 0; i < numSample; ++i)
            availableSample.push_front(i);
        return train(dataset, std::move(availableSample), makeInitialFeatures(dataset.features.getCol()));
    }

    template<Scalar T, DecisionTreeType Type>
    std::forward_list<size_t> DecisionTree<T, Type>::makeInitialFeatures(size_t numFeature) {
        std::forward_list<size_t> availableFeature;
        for (size_t i = 0; i < numFeature; ++i)
            availableFeature.push_front(i);
        return availableFeature;
    }

    template<Scalar T, DecisionTreeType Type>
    T DecisionTree<T, Type>::checkStopRecursion(
            const Dataset& dataset,
            const std::forward_list<size_t>& availableSample,
            const std::forward_list<size_t>& availableFeature) {
        assert(!availableSample.empty());
        constexpr bool isClassifyTree = Type == DecisionTreeType::Classify;
        const auto& features = dataset.features;
        const auto& labels = dataset.labels;

        if constexpr (isClassifyTree) {
            const T label = labels[availableSample.front()];
            bool isAllLabelSame = true;
            for (size_t index : availableSample)
                isAllLabelSame &= label == labels[index];
            
            if (isAllLabelSame)
                return label;
        }

        const bool isFeatureEmpty = availableFeature.empty();
        if (isFeatureEmpty) {
            T prediction;
            if constexpr (isClassifyTree)
                prediction = findCommonLabel(labels, availableSample);
            else
                prediction = makeAverageLabel(labels, availableSample);
            return prediction;
        }

        for (auto featureId : availableFeature) {
            const T feature = features(availableSample.front(), featureId);
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

    template<Scalar T, DecisionTreeType Type>
    DecisionTree<T, Type> DecisionTree<T, Type>::doRecursion(
            const Dataset& dataset,
            const std::forward_list<size_t>& availableSample,
            const std::forward_list<size_t>& availableFeature,
            size_t featureId,
            T splitPoint,
            std::invocable<const Dataset&, std::forward_list<size_t>, std::forward_list<size_t>&> auto fn) {
        std::forward_list<size_t> newAvailableFeature{};
        for (auto feature : availableFeature) {
            if (feature != featureId)
                newAvailableFeature.push_front(feature);
        }

        Array<DecisionTree> subTrees{};
        std::forward_list<size_t> list1{}, list2{};
        subTrees.reserve(2U);
        const bool isContinuous = dataset.isFeatureContinuous[featureId];
        if (isContinuous) {
            for (auto sample : availableSample) {
                if (dataset.features(sample, featureId) < splitPoint)
                    list1.push_front(sample);
                else
                    list2.push_front(sample);
            }
        }
        else {
            for (auto sample : availableSample) {
                if (dataset.features(sample, featureId) == splitPoint)
                    list1.push_front(sample);
                else
                    list2.push_front(sample);
            }
        }
        subTrees.append(fn(dataset, std::move(list1), newAvailableFeature));
        subTrees.append(fn(dataset, std::move(list2), newAvailableFeature));
        return DecisionTree(featureId, splitPoint, std::move(subTrees), isContinuous ? Regression : Classify);
    }

    template<Scalar T, DecisionTreeType Type>
    DecisionTree<T, Type> DecisionTree<T, Type>::train(
            const Dataset& dataset,
            std::forward_list<size_t> availableSample,
            std::forward_list<size_t> availableFeature) {
        using TrainFunctor = DecisionTree (*)(const Dataset& dataset, std::forward_list<size_t>, std::forward_list<size_t>);
        const size_t numFeature = dataset.features.getCol();
        const T criteria = checkStopRecursion(dataset, availableSample, availableFeature);
        const bool shouldStopRecursion = !std::isnan(double(criteria));
        if (shouldStopRecursion)
            return DecisionTree(numFeature, criteria, {}, NodeType::Classify);

        auto pair = selectOptimalFeature(dataset, availableSample, availableFeature, getLossFunctor());
        return doRecursion<TrainFunctor>(dataset, availableSample, availableFeature, pair.first, std::move(pair.second), train);
    }


    template<Scalar T, DecisionTreeType Type>
    T DecisionTree<T, Type>::findCommonLabel(
            const VectorType& labels,
            const std::forward_list<size_t>& availableSample) {
        assert(!availableSample.empty());
        std::multiset<T> set{};
        for (size_t sample : availableSample)
            set.insert(labels[sample]);
        
        size_t numCount = 0;
        T result = 0;
        for (auto ite = set.begin(); ite != set.end(); ++ite) {
            size_t temp = set.count(*ite);
            if (temp > numCount) {
                numCount = temp;
                result = *ite;
            }
        }
        return result;
    }

    template<Scalar T, DecisionTreeType Type>
    T DecisionTree<T, Type>::makeAverageLabel(const VectorType& labels, const std::forward_list<size_t>& availableSample) {
        T result = 0;
        size_t count = 0;
        for (auto sample : availableSample) {
            result += labels[sample];
            count += 1;
        }
        result /= T(count);
        return result;
    }

    template<Scalar T, DecisionTreeType Type>
    T DecisionTree<T, Type>::giniIndex(const Dataset& dataset, const Array<size_t>& availableSample) {
        std::multiset<T> set{};
        for (size_t sample : availableSample)
            set.insert(dataset.labels[sample]);
        
        T result = 1;
        const T numSample = std::distance(availableSample.begin(), availableSample.end());
        for (auto ite = set.begin(); ite != set.end(); ++ite) {
            const T temp = set.count(*ite);
            result -= square(temp / numSample);
        }
        return result;
    }

    template<Scalar T, DecisionTreeType Type>
    T DecisionTree<T, Type>::mse(const Dataset& dataset, const Array<size_t>& availableSample) {
        T sum = 0;
        T sumSquare = 0;
        size_t count = 0;
        for (auto sample : availableSample) {
            const T label = dataset.labels[sample];
            sum += label;
            sumSquare += square(label);
            count += 1;
        }
        const T factor = reciprocal(T(count));
        return sumSquare * factor - square(sum * factor);
    }

    template<Scalar T, DecisionTreeType Type>
    std::pair<size_t, T>
    DecisionTree<T, Type>::selectOptimalFeature(
            const Dataset& dataset,
            const std::forward_list<size_t>& availableSample,
            const std::forward_list<size_t>& availableFeature,
            LossFunctor functor) {
        assert(!availableFeature.empty());

        const size_t numAvailableSample = std::distance(availableSample.begin(), availableSample.end());
        size_t optimalFeatureId = dataset.features.getCol();
        T optimalSplitPoint = 0;
        T minLoss = std::numeric_limits<T>::max();
        Array<size_t> list1{}, list2{};
        list1.reserve(numAvailableSample);
        list2.reserve(numAvailableSample);
        for (auto featureId : availableFeature) {
            VectorType featureVector{};
            {
                std::set<T> featureSet{};
                for (auto sample : availableSample)
                    featureSet.insert(dataset.features(sample, featureId));

                featureVector.resize(featureSet.size());
                size_t i = 0;
                for (auto ite = featureSet.begin(); ite != featureSet.end(); ++ite) {
                    featureVector[i] = *ite;
                    i += 1;
                }
            }

            T splitPoint = 0;
            T loss = std::numeric_limits<T>::max();
            if (dataset.isFeatureContinuous[featureId]) {
                for (size_t i = 0; i < featureVector.getLength() - 1; ++i) {
                    const T midpoint = (featureVector[i] + featureVector[i + 1]) * 0.5;
                    size_t weight1 = 0;
                    for (auto sample : availableSample) {
                        if (dataset.features(sample, featureId) < midpoint) {
                            list1.grow(sample);
                            weight1 += 1;
                        }
                        else
                            list2.grow(sample);
                    }

                    const size_t weight2 = featureVector.getLength() - weight1;
                    const T temp = functor(dataset, list1) * T(weight1) + functor(dataset, list2) * T(weight2);
                    if (temp < loss) {
                        loss = temp;
                        splitPoint = midpoint;
                    }
                    list1.clear();
                    list2.clear();
                }
            }
            else {
                for (size_t i = 0; i < featureVector.getLength(); ++i) {
                    size_t weight1 = 0;
                    for (auto sample : availableSample) {
                        if (dataset.features(sample, featureId) == featureVector[i]) {
                            list1.grow(sample);
                            weight1 += 1;
                        }
                        else
                            list2.grow(sample);
                    }

                    const size_t weight2 = featureVector.getLength() - weight1;
                    const T temp = functor(dataset, list1) * T(weight1) + functor(dataset, list2) * T(weight2);
                    if (temp < loss) {
                        loss = temp;
                        splitPoint = featureVector[i];
                    }
                    list1.clear();
                    list2.clear();
                }
            }

            if (loss <= minLoss) {
                minLoss = loss;
                optimalFeatureId = featureId;
                optimalSplitPoint = splitPoint;
            }
        }
        return {optimalFeatureId, optimalSplitPoint};
    }
}
