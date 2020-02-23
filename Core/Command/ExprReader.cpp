#include <stack>
#include "../Header/ExprReader.h"
#include "../Header/Const.h"

extern const Const_1 const_1;

ExprReader::ExprReader(const std::string& s) {
    std::stack<char> stack{};
    std::string temp;
    for(auto c : s) {
        switch(c) {
            case '-':
                clearString(temp);
                while(!stack.empty() && (stack.top() == '*' || stack.top() == '/')) {
                    std::string temp1;
                    temp1.push_back(stack.top());
                    anti_poland.push_back(temp1);
                    stack.pop();
                }
                if(!anti_poland.empty())
                    stack.push('+');
                temp.push_back(c);
                break;
            case '+':
                clearString(temp);
                while(!stack.empty() && (stack.top() == '*' || stack.top() == '/')) {
                    std::string temp1;
                    temp1.push_back(stack.top());
                    anti_poland.push_back(temp1);
                    stack.pop();
                }
                if(!anti_poland.empty())
                    stack.push('+');
                break;
            case '*':
            case '/':
            case '(':
                clearString(temp);
                stack.push(c);
                break;
            case ')':
                clearString(temp);
                while(!stack.empty()) {
                    if(stack.top() == '(') {
                        stack.pop();
                        break;
                    }
                    std::string temp1;
                    temp1.push_back(stack.top());
                    anti_poland.push_back(temp1);
                    stack.pop();
                }
                break;
            case ' ':
                break;
            default:
                temp.push_back(c);
        }
    }
    anti_poland.push_back(temp);

    while(!stack.empty()) {
        std::string temp1;
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
            case '*':
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
                    stack.push(*n1 / *n2);
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
//Only be used in construct function ExprReader(). Check the temp string and store the useful data.
void ExprReader::clearString(std::string& s) {
    if(!s.empty()) {
        anti_poland.push_back(s);
        s.clear();
    }
}