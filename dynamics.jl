# Walking Arm Trebuchet Optimization Program
# Copyright (C) 2021  Max Bareiss

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

module Falling

using Symbolics
using Latexify

@variables t Im1 IB IL Im2 m1 mB mL m2 L1 L2 L3 β1 β2 θ1(t) θ2(t) θ3(t) h0 g

D = Differential(t)

rL = sqrt((h0/2)^2 + L2^2 + L2*h0*cos(β2))
θL = acos((L2^2+rL^2-(h0/2)^2)/(2*L2*rL))

T = 1/2*m1*((L1+L2)*D(θ1))^2 + 
    1/2*m1*((L1+L2)/2*D(θ1))^2 +
    1/2*mL*(rL*D(θ1))^2 +
    1/2*m2*(L3*D(θ1))^2 +
    D(θ1)^2*1/2*(Im1+IB+IL+Im2)

V = m1*sin(θ1)*(L1+L2)*g + mB*sin(θ1)*((L1+L2)/2)*g + mL*sin(θ1+θL)*rL*g + m2*sin(θ1+β1)*g

L = T-V

#print(latexify(L))

eqn1 = D(Symbolics.derivative(L,D(θ1))) - Symbolics.derivative(L,θ1) == 0

print(latexify(simplify(expand_derivatives(eqn1),expand=true)))

end

module Swinging

using Symbolics
using Latexify

@variables t Im1 IB IL Im2 m1 mB mL m2 L1 L2 L3 β1 β2 θ1(t) θ2(t) θ3(t) h0 g

D = Differential(t)

xB = h0*cos(θ2) + (L1-L2)/2*cos(θ1)
yB = h0*sin(θ2) + (L1-L2)/2*sin(θ1)
vBx = D(xB)
vBy = D(yB)

xm1 = h0*cos(θ2) + L1*cos(θ1)
ym1 = h0*sin(θ2) + L1*sin(θ1)
vm1x = D(xm1)
vm1y = D(ym1)

xm2 = h0*cos(θ2) - L2*cos(θ1) + L3*cos(θ1 + β1)
ym2 = h0*sin(θ2) - L2*sin(θ1) + L3*sin(θ1 + β1)
vm2x = D(xm2)
vm2y = D(ym2)

T = 1/2*mL*(h0/2*D(θ2))^2 +
    1/2*mB*(vBx^2+vBy^2) +
    1/2*m1*(vm1x^2+vm1y^2) +
    1/2*m2*(vm2x^2+vm2y^2) +
    1/2*IL*D(θ2)^2 + 
    1/2*(IB+Im1+Im2)*D(θ1)^2

V = mL*g*h0/2*sin(θ2) + mB*g*yB + m1*g*ym1 + m2*g*ym2

L = T-V

print("\n\nEQUATION 1\n")
eqn1 = D(Symbolics.derivative(L,D(θ1))) - Symbolics.derivative(L,θ1) == 0
print(latexify(simplify(expand_derivatives(eqn1),expand=true,threaded=true)))

print("\n\nEQUATION 2\n")
eqn2 = D(Symbolics.derivative(L,D(θ2))) - Symbolics.derivative(L,θ2) == 0
print(latexify(simplify(expand_derivatives(eqn2),expand=true,threaded=false)))

end

module SwingingBall

using Symbolics
using Latexify

@variables t Im1 IB IL Im2 m1 mB mL m2 L1 L2 L3 β1 β2 θ1(t) θ2(t) θ3(t) h0 g

D = Differential(t)

xB = h0*cos(θ2) + (L1-L2)/2*cos(θ1)
yB = h0*sin(θ2) + (L1-L2)/2*sin(θ1)
vBx = D(xB)
vBy = D(yB)

xm1 = h0*cos(θ2) + L1*cos(θ1)
ym1 = h0*sin(θ2) + L1*sin(θ1)
vm1x = D(xm1)
vm1y = D(ym1)

xm2 = h0*cos(θ2) - L2*cos(θ1) + L3*cos(θ3)
ym2 = h0*sin(θ2) - L2*sin(θ1) + L3*sin(θ3)
vm2x = D(xm2)
vm2y = D(ym2)

T = 1/2*mL*(h0/2*D(θ2))^2 +
    1/2*mB*(vBx^2+vBy^2) +
    1/2*m1*(vm1x^2+vm1y^2) +
    1/2*m2*(vm2x^2+vm2y^2) +
    1/2*IL*D(θ2)^2 + 
    1/2*(IB+Im1)*D(θ1)^2 +
    1/2*Im2*D(θ3)^2

V = mL*g*h0/2*sin(θ2) + mB*g*yB + m1*g*ym1 + m2*g*ym2

L = T-V

print("\n\nEQUATION 1\n")
eqn1 = D(Symbolics.derivative(L,D(θ1))) - Symbolics.derivative(L,θ1) == 0
print(simplify(expand_derivatives(eqn1),expand=true,threaded=false))

print("\n\nEQUATION 2\n")
eqn2 = D(Symbolics.derivative(L,D(θ2))) - Symbolics.derivative(L,θ2) == 0
print(simplify(expand_derivatives(eqn2),expand=true,threaded=false))

print("\n\nEQUATION 3\n")
eqn2 = D(Symbolics.derivative(L,D(θ3))) - Symbolics.derivative(L,θ3) == 0
print(simplify(expand_derivatives(eqn2),expand=true,threaded=false))

print("\n")

end
