#include <iostream>
#include <ctime>
#include "Header/Cube2D.h"
/*
 * This is a program used to solve 2D Cube.
 * It is a easter egg.
 *
 * Warning: Do not run more than 10 steps.
 */
QHash<char*, char*>* temp;

void solveCube() {
    //Do not use auto to avoid crash.
    auto hash = new QHash<char*, char*>;
    temp = new QHash<char*, char*>;
    auto iterator = new QHash<char*, char*>::const_iterator;
    initialize(hash);

    bool goon = true;
    time_t start;
    for(int i = 1; i <= maxStep; i++) {
        start = time(nullptr);
        calculate(hash, iterator);
        std::cout << "The " << i << "th calculate finished, using " << (time(nullptr) - start) << "s. Press 1 to continue." << std::endl;
        std::cin >> goon;
        if(!goon)
            exit(0);
        hash = temp;
        temp = new QHash<char*, char*>;
    }
}

void initialize(QHash<char*, char*>* hash) {
    char* array = new char[length];
    for(int i = 1; i <= length; i++) {
        std::cout << "Input the " << i << "th element." << std::endl;
        std::cin >> array[i - 1];
    }
    hash->insert(new char[maxStep]{0}, array);
}


void calculate(QHash<char*, char*>* hash, QHash<char*, char*>::const_iterator* iterator) {
    for(*iterator = hash->begin(); *iterator != hash->end(); (*iterator)++) {
        char* trans = iterator->key();
        int newTransIndex = 0;
        for(int i = -1; i < maxStep - 1; i++) {
            if(trans[i + 1] == 0) {
                newTransIndex = i + 1;
                break;
            }
        }
        int lastTrans = newTransIndex == 0 ? 0 : trans[newTransIndex - 1];

        if(lastTrans != 4) {
            trans[newTransIndex] = 1;
            transX(iterator);
        }

        if(lastTrans != 5) {
            trans[newTransIndex] = 2;
            transY(iterator);
        }

        if(lastTrans != 6) {
            trans[newTransIndex] = 3;
            transZ(iterator);
        }

        if(lastTrans != 1 && lastTrans != 4) {
            trans[newTransIndex] = 4;
            antiTransX(iterator);
        }

        if(lastTrans != 2 && lastTrans != 5) {
            trans[newTransIndex] = 5;
            antiTransY(iterator);
        }

        if(lastTrans != 3 && lastTrans != 6) {
            trans[newTransIndex] = 6;
            antiTransZ(iterator);
        }
        delete[] iterator->key();
        delete[] iterator->value();
    }
    delete hash;
}
//Trans index: 1
void transX(QHash<char*, char*>::const_iterator* iterator) {
    char* result = new char[length];
    char* trans = new char[maxStep];
    memcpy(result, iterator->value(), length * sizeof(char));
    memcpy(trans, iterator->key(), maxStep * sizeof(char));
    replace(0,1,result);
    replace(1,2,result);
    replace(2,3,result);
    replace(3,0,result);
    replace(10,16,result);
    replace(9,19,result);
    replace(16,14,result);
    replace(19,13,result);
    replace(13,23,result);
    replace(14,20,result);
    replace(20,10,result);
    replace(23,9,result);
    checkAnswer(result, trans);
}
//Trans index: 2
void transY(QHash<char*, char*>::const_iterator* iterator) {
    char* result = new char[length];
    char* trans = new char[maxStep];
    memcpy(result, iterator->value(), length * sizeof(char));
    memcpy(trans, iterator->key(), maxStep * sizeof(char));
    replace(8,9,result);
    replace(9,10,result);
    replace(10,11,result);
    replace(11,8,result);
    replace(7,16,result);
    replace(4,17,result);
    replace(16,3,result);
    replace(17,0,result);
    replace(0,23,result);
    replace(3,22,result);
    replace(22,7,result);
    replace(23,4,result);
    checkAnswer(result, trans);
}
//Trans index: 3
void transZ(QHash<char*, char*>::const_iterator* iterator) {
    char* result = new char[length];
    char* trans = new char[maxStep];
    memcpy(result, iterator->value(), length * sizeof(char));
    memcpy(trans, iterator->key(), maxStep * sizeof(char));
    replace(16,17,result);
    replace(17,18,result);
    replace(18,19,result);
    replace(19,16,result);
    replace(0,8,result);
    replace(1,9,result);
    replace(8,6,result);
    replace(9,7,result);
    replace(6,14,result);
    replace(7,15,result);
    replace(14,0,result);
    replace(15,1,result);
    checkAnswer(result, trans);
}
//Trans index: 4
void antiTransX(QHash<char*, char*>::const_iterator* iterator) {
    char* result = new char[length];
    char* trans = new char[maxStep];
    memcpy(result, iterator->value(), length * sizeof(char));
    memcpy(trans, iterator->key(), maxStep * sizeof(char));
    replace(0,3,result);
    replace(1,0,result);
    replace(2,1,result);
    replace(3,2,result);
    replace(9,23,result);
    replace(10,20,result);
    replace(13,19,result);
    replace(14,16,result);
    replace(16,10,result);
    replace(19,9,result);
    replace(20,14,result);
    replace(23,13,result);
    checkAnswer(result, trans);
}
//Trans index: 5
void antiTransY(QHash<char*, char*>::const_iterator* iterator) {
    char* result = new char[length];
    char* trans = new char[maxStep];
    memcpy(result, iterator->value(), length * sizeof(char));
    memcpy(trans, iterator->key(), maxStep * sizeof(char));
    replace(8,11,result);
    replace(9,8,result);
    replace(10,9,result);
    replace(11,10,result);
    replace(4,23,result);
    replace(7,22,result);
    replace(0,17,result);
    replace(3,16,result);
    replace(17,4,result);
    replace(16,7,result);
    replace(23,0,result);
    replace(22,3,result);
    checkAnswer(result, trans);
}
//Trans index: 6
void antiTransZ(QHash<char*, char*>::const_iterator* iterator) {
    char* result = new char[length];
    char* trans = new char[maxStep];
    memcpy(result, iterator->value(), length * sizeof(char));
    memcpy(trans, iterator->key(), maxStep * sizeof(char));
    replace(16,19,result);
    replace(17,16,result);
    replace(18,17,result);
    replace(19,18,result);
    replace(0,14,result);
    replace(1,15,result);
    replace(7,9,result);
    replace(6,8,result);
    replace(9,1,result);
    replace(8,0,result);
    replace(14,6,result);
    replace(15,7,result);
    checkAnswer(result, trans);
}


bool checkAnswer(char* result, char* trans) {
    bool isAnswer = true;
    for(int i=1; i<5; i++) {
        if(result[4 * i - 4] != result[4 * i - 3] || result[4 * i - 3] != result[4 * i - 2] || result[4 * i - 2] != result[4 * i - 1]) {
            isAnswer = false;
            break;
        }
    }

    if(isAnswer){
        std::cout << "Get answer";
        for(int i=0; i < maxStep; i++)
            std::cout << trans[i];
        std::cout << "\n";
        exit(0);
    }
    else {
        temp->insert(trans, result);
    }
}

void replace(int i,int j,char* matrix) {
    int tempInt = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = char(tempInt);
}