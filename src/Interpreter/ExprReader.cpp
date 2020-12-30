/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
#include <stack>
#include "Physica/Interpreter/ExprReader.h"

using namespace Physica::Interpreter;

ExprReader::ExprReader(const std::wstring& str) {
    size_t pointer = 0;
    std::stack<wchar_t> stack{};
    std::wstring temp;
    while(pointer != str.size()) {
        if(isSign(str[pointer])) {
            bool got_sign = false;
            do {
                if(!got_sign) {
                    switch(str[pointer]) {
                        case '-':
                            temp.push_back('-');
                        case '+':
                            got_sign = true;
                            while(!stack.empty() && (stack.top() == L'×' || stack.top() == '/')) {
                                std::wstring temp1;
                                temp1.push_back(stack.top());
                                anti_poland.push_back(temp1);
                                stack.pop();
                            }
                            if(!anti_poland.empty())
                                stack.push('+');
                            break;
                        case L'×':
                        case '/':
                            got_sign = true;
                            if(!stack.empty() && (stack.top() == L'×' || stack.top() == '/')) {
                                std::wstring temp1{stack.top()};
                                anti_poland.push_back(temp1);
                                stack.pop();
                                stack.push(str[pointer]);
                                break;
                            }
                        case '(':
                            stack.push(str[pointer]);
                            break;
                        case ')':
                            got_sign = true;
                            while(!stack.empty()) {
                                if(stack.top() == '(') {
                                    stack.pop();
                                    break;
                                }
                                std::wstring temp1;
                                temp1.push_back(stack.top());
                                anti_poland.push_back(temp1);
                                stack.pop();
                            }
                            break;
                        default:;
                    }
                }
                ++pointer;
            } while(isSign(str[pointer]));
            //If got_sign is false, we mean we only have ' ' or '(', where we should do multiply.
            if(!got_sign) {
                if(stack.top() == '(') {
                    stack.pop();
                    stack.push(L'×');
                    stack.push('(');
                }
            }
        }
        else {
            while(pointer != str.size() && !isSign(str[pointer])) {
                temp.push_back(str[pointer]);
                ++pointer;
            }
            anti_poland.push_back(temp);
            temp.clear();
        }
    }

    while(!stack.empty()) {
        std::wstring temp1;
        temp1.push_back(stack.top());
        anti_poland.push_back(temp1);
        stack.pop();
    }
}

MultiScalar ExprReader::calc() {
    std::stack<MultiScalar> stack{};
    for(auto& str : anti_poland) {
        switch(str.front()) {
            case '+':
                if(stack.size() > 1) {
                    MultiScalar n1 = stack.top();
                    stack.pop();
                    MultiScalar n2 = stack.top();
                    stack.pop();
                    stack.push(n1 + n2);
                }
                break;
            case L'×':
                if(stack.size() > 1) {
                    MultiScalar n1 = stack.top();
                    stack.pop();
                    MultiScalar n2 = stack.top();
                    stack.pop();
                    stack.push(n1 * n2);
                }
                break;
            case '/':
                if(stack.size() > 1) {
                    MultiScalar n1 = stack.top();
                    stack.pop();
                    MultiScalar n2 = stack.top();
                    stack.pop();
                    stack.push(n2 / n1);
                }
                break;
            default:
                stack.push(MultiScalar(str.c_str()));
        }
    }
    return stack.top();
}

bool ExprReader::isSign(wchar_t c) {
    switch(c) {
        case '+':
        case '-':
        case L'×':
        case '/':
        case '(':
        case ')':
        case ' ':
            return true;
        default:
            return false;
    }
}