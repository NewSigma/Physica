#ifndef PHYSICA_C_CUBE2_H
#define PHYSICA_C_CUBE2_H

#include <QtCore/QHash>

const int length = 24;
const int maxStep = 14;

void solveCube();
void initialize(QHash<char*, char*>* hash);
void calculate(QHash<char*, char*>* hash, QHash<char*, char*>::const_iterator* iterator);
void replace(int i,int j,char* matrix);
bool checkAnswer(char* result, char* trans);
void transX(QHash<char*, char*>::const_iterator* iterator);
void transY(QHash<char*, char*>::const_iterator* iterator);
void transZ(QHash<char*, char*>::const_iterator* iterator);
void antiTransX(QHash<char*, char*>::const_iterator* iterator);
void antiTransY(QHash<char*, char*>::const_iterator* iterator);
void antiTransZ(QHash<char*, char*>::const_iterator* iterator);

#endif
