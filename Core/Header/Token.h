/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKEN_H
#define PHYSICA_TOKEN_H

class Token {
    const char* data;
public:
    enum TokenType {
        Null,
        Number,
        Identifier,
        Keyword
    } const type;

    Token(const char* data, TokenType type);
    ~Token();
};

#endif
