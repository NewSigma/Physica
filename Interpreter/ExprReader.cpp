#include <stack>
#include "Interpreter/Header/ExprReader.h"
#include "Core/Header/Scalar.h"

using namespace Physica::Interpreter;

ExprReader::ExprReader(const std::wstring& str) {
    int pointer = 0;
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

Scalar ExprReader::calc() {
    std::stack<Scalar> stack{};
    for(auto& str : anti_poland) {
        switch(str.front()) {
            case '+':
                if(stack.size() > 1) {
                    Scalar n1 = stack.top();
                    stack.pop();
                    Scalar n2 = stack.top();
                    stack.pop();
                    stack.push(n1 + n2);
                }
                break;
            case L'×':
                if(stack.size() > 1) {
                    Scalar n1 = stack.top();
                    stack.pop();
                    Scalar n2 = stack.top();
                    stack.pop();
                    stack.push(n1 * n2);
                }
                break;
            case '/':
                if(stack.size() > 1) {
                    Scalar n1 = stack.top();
                    stack.pop();
                    Scalar n2 = stack.top();
                    stack.pop();
                    stack.push(n2 / n1);
                }
                break;
            default:
                stack.push(Scalar(str));
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