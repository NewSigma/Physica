#include <stack>
#include "../Header/ExprReader.h"
#include "../Header/Const.h"

extern const Const_1 const_1;

ExprReader::ExprReader(const std::wstring& str) {
    int pointer = 0;
    std::stack<wchar_t> stack{};
    std::wstring temp;
    while(pointer != str.size()) {
        if(isSign(str[pointer])) {
            bool got_sign = false;
            while(isSign(str[pointer])) {
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
            }
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

RealNumber* ExprReader::calc() {
    std::stack<RealNumber*> stack{};
    for(auto& str : anti_poland) {
        RealNumber* n1 = nullptr, *n2 = nullptr;
        switch(str.front()) {
            case '+':
                if(stack.size() > 1) {
                    n1 = stack.top();
                    stack.pop();
                    n2 = stack.top();
                    stack.pop();
                    stack.push(*n1 + *n2);
                }
                break;
            case L'×':
                if(stack.size() > 1) {
                    n1 = stack.top();
                    stack.pop();
                    n2 = stack.top();
                    stack.pop();
                    stack.push(*n1 * *n2);
                }
                break;
            case '/':
                if(stack.size() > 1) {
                    n1 = stack.top();
                    stack.pop();
                    n2 = stack.top();
                    stack.pop();
                    stack.push(*n2 / *n1);
                }
                break;
            default:
                stack.push(new RealNumber(str));
        }
        delete n1;
        delete n2;
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