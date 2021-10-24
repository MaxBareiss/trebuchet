


if __name__ == "__main__":
    eqns = ['((IB + Im1)*Differential(t)(Differential(t)(θ1(t))) + L1*g*m1*cos(θ1(t)) + m1*(L1^2)*(cos(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) + m1*(L1^2)*(sin(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) + m2*(L2^2)*(cos(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) + m2*(L2^2)*(sin(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) + mB*(h0^2)*(sin(θ2(t))^2)*Differential(t)(Differential(t)(θ2(t))) + 0.5L1*g*mB*cos(θ1(t)) + 0.25mB*(L2^2)*(cos(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) + mB*(h0^2)*(Differential(t)(θ2(t))^2)*cos(θ2(t))*sin(θ2(t)) + 0.25mB*(L1^2)*(cos(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) + L2*L3*m2*(Differential(t)(θ3(t))^2)*cos(θ1(t))*sin(θ3(t)) + L1*h0*m1*(Differential(t)(θ2(t))^2)*sin(θ1(t))*cos(θ2(t)) + L1*h0*m1*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ2(t))) + L1*h0*m1*sin(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ2(t))) + L2*h0*m2*(Differential(t)(θ2(t))^2)*cos(θ1(t))*sin(θ2(t)) + 0.5L1*L2*mB*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ1(t)) + 0.5L1*h0*mB*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ2(t)) + 0.5L1*h0*mB*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ2(t))) + 0.5L1*h0*mB*sin(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ1(t))) + 0.5L2*h0*mB*(Differential(t)(θ2(t))^2)*cos(θ1(t))*sin(θ2(t)) - 0.5L2*g*mB*cos(θ1(t)) - L2*g*m2*cos(θ1(t)) - 0.25mB*(L1^2)*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ1(t)) - 0.25mB*(L2^2)*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ1(t)) - 0.5L1*L2*mB*(cos(θ1(t))^2)*Differential(t)(Differential(t)(θ1(t))) - L2*L3*m2*(Differential(t)(θ3(t))^2)*sin(θ1(t))*cos(θ3(t)) - L2*L3*m2*cos(θ1(t))*cos(θ3(t))*Differential(t)(Differential(t)(θ3(t))) - L1*h0*m1*(Differential(t)(θ2(t))^2)*cos(θ1(t))*sin(θ2(t)) - L2*L3*m2*sin(θ1(t))*sin(θ3(t))*Differential(t)(Differential(t)(θ3(t))) - 0.5L2*h0*mB*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ2(t)) - 0.5L2*h0*mB*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ2(t))) - L2*h0*m2*(Differential(t)(θ2(t))^2)*sin(θ1(t))*cos(θ2(t)) - L2*h0*m2*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ2(t))) - L2*h0*m2*sin(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ2(t))) - 0.5L1*h0*mB*(Differential(t)(θ2(t))^2)*cos(θ1(t))*sin(θ2(t)) - 0.5L2*h0*mB*sin(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ1(t)))) == 0',
            '(IL*Differential(t)(Differential(t)(θ2(t))) + g*h0*m1*cos(θ2(t)) + g*h0*m2*cos(θ2(t)) + g*h0*mB*cos(θ2(t)) + 0.25mL*(h0^2)*Differential(t)(Differential(t)(θ2(t))) + m1*(h0^2)*(cos(θ2(t))^2)*Differential(t)(Differential(t)(θ2(t))) + (h0^2)*(m1 + m2)*(sin(θ2(t))^2)*Differential(t)(Differential(t)(θ2(t))) + (h0^2)*(m2 + mB)*(cos(θ2(t))^2)*Differential(t)(Differential(t)(θ2(t))) + 0.25mB*(L1^2)*(Differential(t)(θ1(t))^3)*cos(θ1(t)) + 0.25mB*(L2^2)*(Differential(t)(θ1(t))^3)*cos(θ1(t)) + 0.5g*h0*mL*cos(θ2(t)) + L2*h0*m2*(Differential(t)(θ1(t))^2)*sin(θ1(t))*cos(θ2(t)) + L1*h0*m1*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ1(t))) + L1*h0*m1*sin(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ1(t))) + L3*h0*m2*cos(θ2(t))*cos(θ3(t))*Differential(t)(Differential(t)(θ3(t))) + 0.5mB*(L1^2)*Differential(t)(θ1(t))*sin(θ1(t))*Differential(t)(Differential(t)(θ1(t))) + L1*h0*m1*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ2(t)) + L3*h0*m2*(Differential(t)(θ3(t))^2)*sin(θ2(t))*cos(θ3(t)) + L3*h0*m2*sin(θ2(t))*sin(θ3(t))*Differential(t)(Differential(t)(θ3(t))) + 0.5mB*(L2^2)*Differential(t)(θ1(t))*sin(θ1(t))*Differential(t)(Differential(t)(θ1(t))) + 0.5L1*L2*mB*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ1(t)) + 0.5L1*h0*mB*(Differential(t)(θ2(t))^2)*Differential(t)(θ1(t))*cos(θ2(t)) + 0.5L1*h0*mB*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ1(t))) + 0.5L1*h0*mB*Differential(t)(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ2(t))) + 0.5L2*h0*mB*(Differential(t)(θ1(t))^2)*sin(θ1(t))*cos(θ2(t)) + 0.5L1*h0*mB*Differential(t)(θ2(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ1(t))) + 0.5L2*h0*mB*Differential(t)(θ1(t))*cos(θ1(t))*Differential(t)(θ2(t))*sin(θ2(t)) - 0.25mB*(L1^2)*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ1(t)) - 0.5L1*L2*mB*(Differential(t)(θ1(t))^3)*cos(θ1(t)) - 0.25mB*(L2^2)*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ1(t)) - mB*(h0^2)*(Differential(t)(θ2(t))^2)*cos(θ2(t))*sin(θ2(t)) - 0.5L1*h0*mB*(Differential(t)(θ1(t))^2)*sin(θ1(t))*cos(θ2(t)) - L2*h0*m2*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ2(t)) - L2*h0*m2*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ1(t))) - L3*h0*m2*(Differential(t)(θ3(t))^2)*cos(θ2(t))*sin(θ3(t)) - 0.5L2*h0*mB*(Differential(t)(θ2(t))^2)*Differential(t)(θ1(t))*cos(θ2(t)) - 0.5L2*h0*mB*Differential(t)(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ2(t))) - 0.5L2*h0*mB*cos(θ1(t))*cos(θ2(t))*Differential(t)(Differential(t)(θ1(t))) - 0.5L2*h0*mB*Differential(t)(θ2(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ1(t))) - L1*L2*mB*Differential(t)(θ1(t))*sin(θ1(t))*Differential(t)(Differential(t)(θ1(t))) - L2*h0*m2*sin(θ1(t))*sin(θ2(t))*Differential(t)(Differential(t)(θ1(t))) - L1*h0*m1*(Differential(t)(θ1(t))^2)*sin(θ1(t))*cos(θ2(t)) - 0.5L1*h0*mB*Differential(t)(θ1(t))*cos(θ1(t))*Differential(t)(θ2(t))*sin(θ2(t))) == 0',
            '(Im2*Differential(t)(Differential(t)(θ3(t))) + L3*g*m2*cos(θ3(t)) + m2*(L3^2)*(cos(θ3(t))^2)*Differential(t)(Differential(t)(θ3(t))) + m2*(L3^2)*(sin(θ3(t))^2)*Differential(t)(Differential(t)(θ3(t))) + L2*L3*m2*(Differential(t)(θ1(t))^2)*sin(θ1(t))*cos(θ3(t)) + L3*h0*m2*(Differential(t)(θ2(t))^2)*cos(θ2(t))*sin(θ3(t)) + L3*h0*m2*cos(θ2(t))*cos(θ3(t))*Differential(t)(Differential(t)(θ2(t))) + L3*h0*m2*sin(θ2(t))*sin(θ3(t))*Differential(t)(Differential(t)(θ2(t))) - L2*L3*m2*(Differential(t)(θ1(t))^2)*cos(θ1(t))*sin(θ3(t)) - L2*L3*m2*cos(θ1(t))*cos(θ3(t))*Differential(t)(Differential(t)(θ1(t))) - L3*h0*m2*(Differential(t)(θ2(t))^2)*sin(θ2(t))*cos(θ3(t)) - L2*L3*m2*sin(θ1(t))*sin(θ3(t))*Differential(t)(Differential(t)(θ1(t)))) == 0']
    for i, eqn in enumerate(eqns):
        eqn = eqn.replace('Differential(t)(Differential(t)(θ1(t)))','DDOT1')
        eqn = eqn.replace('Differential(t)(Differential(t)(θ2(t)))','DDOT2')
        eqn = eqn.replace('Differential(t)(Differential(t)(θ3(t)))','DDOT3')
        eqn = eqn.replace('Differential(t)(θ1(t))','ω1')
        eqn = eqn.replace('Differential(t)(θ2(t))','ω2')
        eqn = eqn.replace('Differential(t)(θ3(t))','ω3')
        eqn = eqn.replace('θ1(t)','θ1')
        eqn = eqn.replace('θ2(t)','θ2')
        eqn = eqn.replace('θ3(t)','θ3')
        eqn = eqn.replace('5','5*')
        #print(eqn)
        parens = 0
        last_char = ''
        sub_exprs = ['+ ']
        for character in eqn:
            if character == "(":
                parens += 1
                if parens > 1:
                    sub_exprs[-1] += character
            elif character == ")":
                if parens > 1:
                    sub_exprs[-1] += character
                parens -= 1
            elif last_char == '=' and character == '=':
                sub_exprs[-1] = sub_exprs[-1].strip()
                break
            elif character == '=':
                pass
            elif character == "+" and parens == 1:
                sub_exprs[-1] = sub_exprs[-1].strip()
                sub_exprs.append("+ ")
            elif character == "-" and parens == 1:
                sub_exprs[-1] = sub_exprs[-1].strip()
                sub_exprs.append("- ")
            elif character == " " and parens == 1:
                pass
            else:
                sub_exprs[-1] += character
            last_char = character

        #print(sub_exprs)

        print("\n\nEQUATION {}\n\n".format(i+1))

        
        th1 = list()
        th2 = list()
        th3 = list()
        other = list()
        for expr in sub_exprs:
            if "DDOT1" in expr:
                expr = expr.replace("*DDOT1","")
                expr = expr.replace("DDOT1*","")
                th1.append(expr)
            elif "DDOT2" in expr:
                expr = expr.replace("*DDOT2","")
                expr = expr.replace("DDOT2*","")
                th2.append(expr)
            elif "DDOT3" in expr:
                expr = expr.replace("*DDOT3","")
                expr = expr.replace("DDOT3*","")
                th3.append(expr)
            else:
                expr = list(expr)
                if expr[0] == "+":
                    expr[0] = "-"
                elif expr[0] == "-":
                    expr[0] = "+"
                expr = "".join(expr)
                other.append(expr)

        print("DDOT1")
        print(" ".join(th1))

        print("DDOT2")
        print(" ".join(th2))

        print("DDOT3")
        print(" ".join(th3))

        print("OTHER")
        print(" ".join(other))