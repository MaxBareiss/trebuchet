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

module BeamTrebuchet

using DifferentialEquations
using Unitful

Unitful.preferunits(u"ft")

m1 = ustrip(u"kg",10.0u"kg")
mL = ustrip(u"kg",1.0u"kg")
mB = ustrip(u"kg",5.0u"kg")
m2 = ustrip(u"kg",0.2u"kg")
h0 = ustrip(u"m",20.5u"inch")
L1 = ustrip(u"m",14.75u"inch")
L2 = ustrip(u"m",20.0u"inch")
L3 = ustrip(u"m",20.0u"inch")
β1 = ustrip(u"rad",20.0u"°")
β2 = ustrip(u"rad",50.0u"°")
β3 = ustrip(u"rad",0.0u"°")
β4 = ustrip(u"rad",-40.0u"°")
IB = 1/12*mB*(L1+L2)^2 # ustrip(u"kg*m^2",
IL = 1/12*mL*h0^2
Im1 = 2/5*m1*0.2^2
Im2 = 2/5*m2*0.2^2
g = ustrip(u"m/s^2",9.81u"m/s^2")
rL = sqrt((h0/2)^2 + L2^2 + L2*h0*cos(β2))
θL = acos((L2^2+rL^2-(h0/2)^2)/(2*L2*rL))

function falling!(du,u,p,t)
    # u is [θ1, ω1]
    θ1,ω1 = u

    du[1] = u[2]
    du[2] = (-m2*g*cos(θ1+β1)*L3 - g*m1*(L1+L2)*cos(θ1) - g*mL*cos(θ1+θL)*rL - (L1+L2)/2*mB*cos(θ1)*g)/
            (IB + IL + Im1 + Im2 + L3^2*m2 + rL^2*mL + 5/4*L1^2*m1 + 5/4*L2^2*m1 + 5/2*L2*L1*m1)
    evt = θ1 - atan((h0*sin(β2))/(L2+h0*cos(β2))) - ustrip(u"rad",90u"°")
    r = [θ1,ω1,θ1+β2,ω1,θ1+β1,ω1]
end

function falling_cond(u, t, integrator)
    θ1,ω1 = u
    #evt = θ1 - atan((h0*sin(β2))/(L2+h0*cos(β2))) - ustrip(u"rad",90u"°")
    evt = sin(θ1)*L2+sin(θ1+β2)*h0
end

function terminator!(integrator)
    terminate!(integrator)
end

u0 = [ustrip(u"rad",90u"°"),ustrip(u"rad/s",10u"°/s")]

fallingCB = ContinuousCallback(falling_cond,terminator!)

fallingProb = ODEProblem(falling!,u0,[0,10.0])

fallingsol = solve(fallingProb; callback=fallingCB)

# State variable is θ1 ω1 θ2 ω2
u0 = [fallingsol.u[end][1], fallingsol.u[end][2], fallingsol.u[end][1] + β2, fallingsol.u[end][2]]

rL = cos(fallingsol.u[end][1])*L2+cos(fallingsol.u[end][1]+β2)*h0

function swinging!(du,u,p,t)
    # u is θ1 ω1 θ2 ω2
    θ1,ω1,θ2,ω2 = u

    du[1] = u[2]
    du[2] = 5
    du[3] = u[4]
    du[4] = 5

    M = zeros(2,2)

    Q = zeros(2,1)

    M[1,1] = IB + Im1 + Im2 + L1^2*m1 + L3^2*m2 + L2^2*m2 + 1/2*mB*cos(θ1)^2*L2^2 + 1/2*mB*cos(θ1)^2*L1^2 -
             1/4*L1^2*mB*cos(θ1)*sin(θ1) - 1/2*cos(θ1)^2*L1*L2*mB - 2*L2*L3*m2*sin(θ1 + β1)*sin(θ1) -
             2*L2*L3*cos(θ1)*cos(θ1+β1) + 1/2*L1*h0*mB*sin(θ1)*sin(θ2) - 1/2*L2*h0*mB*sin(θ1)*sin(θ2)

    M[1,2] = sin(θ2)^2*h0^2*mB + L3*h0*m2*sin(θ1+β1)*sin(θ2) + L1*h0*m1*cos(θ1)*cos(θ2) + L1*h0*mL*sin(θ1)*sin(θ2) +
             1/2*L2*h0*mB*cos(θ1)*cos(θ2) - L2*h0*m2*cos(θ1)*cos(θ2) - L2*h0*m2*sin(θ1)*sin(θ2) +1/2*L1*h0*mB*cos(θ1)*cos(θ2)

    Q[1,1] = -cos(θ1)*(L1*g*m1+1/2*L1*g*mB-1/2*L2*g*mB-L2*m2) - cos(θ1+β1)*(L3*g*m2) - ω2^2*h0^2*mB*cos(θ2)*sin(θ2) +
             1/4*ω1^2*L2^2*mB*cos(θ1)*sin(θ1) - ω2^2*L1*h0*m1*sin(θ1)*cos(θ2) - ω2^2*L3*h0*m2*sin(θ1+β1)*cos(θ2) - ω2*L3*h0*m2*cos(θ1+β1)*cos(θ2) -
             ω2^2*L2*h0*m2*cos(θ1)*sin(θ2) + ω2^2*L3*h0*m2*cos(θ1+β1)*sin(θ2) + ω2^2*L1*h0*m1*cos(θ1)*sin(θ2) - 1/2*ω1^2*L1*L2*mB*cos(θ1)*sin(θ1) -
             1/2*ω1^2*L1*h0*mB*cos(θ1)*sin(θ2) + 1/2*ω1^2*L2*h0*mB*cos(θ1)*sin(θ2) + ω2^2*L2*h0*m2*sin(θ1)*cos(θ2) + 1/2*ω2^2*L1*h0*mB*cos(θ1)*sin(θ2) -
             1/2*ω2^2*L2*h0*mB*cos(θ1)*sin(θ2)

    M[2,1] = L1*m1*h0*cos(θ1)*cos(θ2) + L1*h0*m1*sin(θ1)*sin(θ2) + 1/2*L1^2*mB*ω1*sin(θ1) + L3*h0*m2*cos(θ1+β1)*cos(θ2) +
             1/2*L2^2*mB*ω1*sin(θ1) + L3*h0*m2*sin(θ1+β1)*sin(θ2) - 1/2*L2*h0*mB*cos(θ1)*cos(θ2) - 1/2*L2*h0*mB*ω2*sin(θ2) -
             L2*h0*m2*cos(θ1)*cos(θ2) - L1*L2*mB*ω1*sin(θ1) - L2*h0*m2*sin(θ1)*sin(θ2) + 1/2*L1*h0*mB*cos(θ1)*cos(θ2) +
             1/2*L1*h0*mB*ω2*sin(θ2)

    M[2,2] = IL + 1/4*h0^2*mL + cos(θ2)^2*h0^2*m1 + sin(θ2)^2*h0^2*(m1+m2) +cos(θ2)^2*h0^2*(m2+mB) - 1/2*L2*h0*mB*ω1*sin(θ2) +
             1/2*L1*h0*mB*ω1*sin(θ2)

    Q[2,1] = -g*h0*cos(θ2)*(m1+m2+mB+1/2*mL) - 1/4*ω1^3*mB*cos(θ1)*(L1^2 + L2^2) + 1/4*ω1^2*L1^2*mB*cos(θ1)*sin(θ1) + 1/2*ω1^3*L1*L2*mB*cos(θ1) -
             ω1^2*L2*h0*m2*sin(θ1)*cos(θ2) + 1/4*ω1^2*L2^2*mB*cos(θ1)*sin(θ1) + ω2^2*h0^2*mB*cos(θ2)*sin(θ2) - ω1^2*L3*h0*m2*cos(θ1+β1)*sin(θ2) -
             θ1^2*L1*h0*m1*cos(θ1)*sin(θ2) - 1/2*ω1^2*L1*L2*mB*cos(θ1)*sin(θ1) + 1/2*ω1^2*L1*h0*mB*sin(θ1)*cos(θ2) + 1/2*ω2^2*L2*h0*mB*ω1*cos(θ2) +
             ω1^2*L3*h0*m2*sin(θ1+β1)*cos(θ2) + ω1^2*L2*h0*m2*cos(θ1)*sin(θ2) - 1/2*ω2^2*L1*h0*mB*ω1*cos(θ2) - 1/2*ω1^2*L2*h0*mB*sin(θ1)*cos(θ2) +
             ω1^2*L1*h0*m1*sin(θ1)*cos(θ2) + 1/2*L1*h0*mB*ω1*ω2*cos(θ1)*sin(θ2) - 1/2*L2*h0*mB*ω1*ω2*cos(θ1)*sin(θ2)

    p = M\Q

    du[2] = p[1]
    du[4] = p[2]

    #evt = θ1 - atan((h0*sin(β2))/(L2+h0*cos(β2))) - ustrip(u"rad",90u"°")
    r = [θ1,ω1,θ2,ω2,θ1+β1,ω1]
end

function swinging_cond(u, t, integrator)
    # u is θ1 ω1 θ2 ω2
    θ1,ω1,θ2,ω2 = u
    evt = θ1 - π + β3
end

swingingCB = ContinuousCallback(swinging_cond,terminator!)

swingingProb = ODEProblem(swinging!,u0,[fallingsol.t[end],30])

swingingsol = solve(swingingProb, Tsit5(), reltol=1e-8, abstol=1e-8; callback=swingingCB)

#print(swingingsol)

function swingingBall!(du,u,p,t)
    # u is θ1 ω1 θ2 ω2 θ3 ω3
    θ1,ω1,θ2,ω2,θ3,ω3 = u

    du[1] = u[2]
    du[2] = 5
    du[3] = u[4]
    du[4] = 5
    du[5] = u[6]
    du[6] = 5

    M = zeros(3,3)

    Q = zeros(3,1)

    M[1,1] = (IB + Im1) + m1*(L1^2)*(cos(θ1)^2) + m1*(L1^2)*(sin(θ1)^2) + m2*(L2^2)*(cos(θ1)^2) + m2*(L2^2)*(sin(θ1)^2) + 
             0.25*mB*(L2^2)*(cos(θ1)^2) + 0.25*mB*(L1^2)*(cos(θ1)^2) + 0.5*L1*h0*mB*sin(θ1)*sin(θ2) - 0.5*L1*L2*mB*(cos(θ1)^2) - 
             0.5*L2*h0*mB*sin(θ1)*sin(θ2)

    M[1,2] = mB*(h0^2)*(sin(θ2)^2) + L1*h0*m1*cos(θ1)*cos(θ2) + L1*h0*m1*sin(θ1)*sin(θ2) + 0.5*L1*h0*mB*cos(θ1)*cos(θ2) -
             0.5*L2*h0*mB*cos(θ1)*cos(θ2) - L2*h0*m2*cos(θ1)*cos(θ2) - L2*h0*m2*sin(θ1)*sin(θ2)

    M[1,3] = L2*L3*m2*cos(θ1)*cos(θ3) - L2*L3*m2*sin(θ1)*sin(θ3)

    Q[1,1] = - L1*g*m1*cos(θ1) - 0.5*L1*g*mB*cos(θ1) - mB*(h0^2)*(ω2^2)*cos(θ2)*sin(θ2) - L2*L3*m2*(ω3^2)*cos(θ1)*sin(θ3) - 
               L1*h0*m1*(ω2^2)*sin(θ1)*cos(θ2) - L2*h0*m2*(ω2^2)*cos(θ1)*sin(θ2) - 0.5*L1*L2*mB*(ω1^2)*cos(θ1)*sin(θ1) -
               0.5*L1*h0*mB*(ω1^2)*cos(θ1)*sin(θ2) - 0.5*L2*h0*mB*(ω2^2)*cos(θ1)*sin(θ2) + 0.5*L2*g*mB*cos(θ1) + L2*g*m2*cos(θ1) +
               0.25*mB*(L1^2)*(ω1^2)*cos(θ1)*sin(θ1) + 0.25*mB*(L2^2)*(ω1^2)*cos(θ1)*sin(θ1) + L2*L3*m2*(ω3^2)*sin(θ1)*cos(θ3) +
               L1*h0*m1*(ω2^2)*cos(θ1)*sin(θ2) + 0.5*L2*h0*mB*(ω1^2)*cos(θ1)*sin(θ2) + L2*h0*m2*(ω2^2)*sin(θ1)*cos(θ2) +
               0.5*L1*h0*mB*(ω2^2)*cos(θ1)*sin(θ2)

    M[2,1] = L1*h0*m1*cos(θ1)*cos(θ2) + L1*h0*m1*sin(θ1)*sin(θ2) + 0.5*mB*(L1^2)*ω1*sin(θ1) + 0.5*mB*(L2^2)*ω1*sin(θ1) +
             0.5*L1*h0*mB*cos(θ1)*cos(θ2) + 0.5*L1*h0*mB*ω2*sin(θ2) - L2*h0*m2*cos(θ1)*cos(θ2) - 0.5*L2*h0*mB*cos(θ1)*cos(θ2) -
             0.5*L2*h0*mB*ω2*sin(θ2) - L1*L2*mB*ω1*sin(θ1) - L2*h0*m2*sin(θ1)*sin(θ2)

    M[2,2] = IL + 0.25*mL*(h0^2) + m1*(h0^2)*(cos(θ2)^2) + (h0^2)*(m1 + m2)*(sin(θ2)^2) + (h0^2)*(m2 + mB)*(cos(θ2)^2) +
             0.5*L1*h0*mB*ω1*sin(θ2) - 0.5*L2*h0*mB*ω1*sin(θ2)

    M[2,3] = L3*h0*m2*cos(θ2)*cos(θ3) + L3*h0*m2*sin(θ2)*sin(θ3)

    Q[2,1] = - g*h0*m1*cos(θ2) - g*h0*m2*cos(θ2) - g*h0*mB*cos(θ2) - 0.25*mB*(L1^2)*(ω1^3)*cos(θ1) - 0.25*mB*(L2^2)*(ω1^3)*cos(θ1) -
               0.5*g*h0*mL*cos(θ2) - L2*h0*m2*(ω1^2)*sin(θ1)*cos(θ2) - L1*h0*m1*(ω1^2)*cos(θ1)*sin(θ2) - L3*h0*m2*(ω3^2)*sin(θ2)*cos(θ3) -
               0.5*L1*L2*mB*(ω1^2)*cos(θ1)*sin(θ1) - 0.5*L1*h0*mB*(ω2^2)*ω1*cos(θ2) - 0.5*L2*h0*mB*(ω1^2)*sin(θ1)*cos(θ2) -
               0.5*L2*h0*mB*ω1*cos(θ1)*ω2*sin(θ2) + 0.25*mB*(L1^2)*(ω1^2)*cos(θ1)*sin(θ1) + 0.5*L1*L2*mB*(ω1^3)*cos(θ1) +
               0.25*mB*(L2^2)*(ω1^2)*cos(θ1)*sin(θ1) + mB*(h0^2)*(ω2^2)*cos(θ2)*sin(θ2) + 0.5*L1*h0*mB*(ω1^2)*sin(θ1)*cos(θ2) +
               L2*h0*m2*(ω1^2)*cos(θ1)*sin(θ2) + L3*h0*m2*(ω3^2)*cos(θ2)*sin(θ3) + 0.5*L2*h0*mB*(ω2^2)*ω1*cos(θ2) +
               L1*h0*m1*(ω1^2)*sin(θ1)*cos(θ2) + 0.5*L1*h0*mB*ω1*cos(θ1)*ω2*sin(θ2)

    M[3,1] = - L2*L3*m2*cos(θ1)*cos(θ3) - L2*L3*m2*sin(θ1)*sin(θ3)

    M[3,2] = L3*h0*m2*cos(θ2)*cos(θ3) + L3*h0*m2*sin(θ2)*sin(θ3)

    M[3,3] = Im2 + m2*(L3^2)*(cos(θ3)^2) + m2*(L3^2)*(sin(θ3)^2)

    Q[3,1] = - L3*g*m2*cos(θ3) - L2*L3*m2*(ω1^2)*sin(θ1)*cos(θ3) - L3*h0*m2*(ω2^2)*cos(θ2)*sin(θ3) + L2*L3*m2*(ω1^2)*cos(θ1)*sin(θ3) + 
               L3*h0*m2*(ω2^2)*sin(θ2)*cos(θ3)


    p = M\Q

    du[2] = p[1]
    du[4] = p[2]
    du[6] = p[3]

    r = [θ1,ω1,θ2,ω2,θ3,ω3]
end

function launch_cond(u, t, integrator)
    # u is θ1 ω1 θ2 ω2
    θ1,ω1,θ2,ω2,θ3,ω3 = u
    evt = θ1 + π + - θ3 + β4
end

# State variable is θ1 ω1 θ2 ω2
u0 = [swingingsol.u[end][1], swingingsol.u[end][2], swingingsol.u[end][3], swingingsol.u[end][4],swingingsol.u[end][1] + β1,swingingsol.u[end][2]]

swingingBallCB = ContinuousCallback(launch_cond,terminator!)

swingingBallProb = ODEProblem(swingingBall!,u0,[swingingsol.t[end],30])

swingingBallsol = solve(swingingBallProb, Tsit5(), reltol=1e-6, abstol=1e-6, dtmin=0.0001, force_dtmin=true; callback=swingingBallCB)

#print(swingingBallsol)

step1_start = fallingsol.t[1]
step1_end = fallingsol.t[end]
step2_start = step1_end
step2_end = swingingsol.t[end]
step3_start = step2_end
step3_end = swingingBallsol.t[end]

t = Vector(range(0,stop=step3_end,length=1000))

step1 = [[x[1],x[2],x[1]+β2,x[2],x[1]+β1,x[2]] for x in fallingsol(t[(t .>= step1_start) .& (t .< step1_end)])]
step2 = [[x[1],x[2],x[3]   ,x[4],x[1]+β1,x[2]] for x in swingingsol(t[(t .>= step2_start) .& (t .< step2_end)])]
step3 = [[x[1],x[2],x[3]   ,x[4],x[5]   ,x[6]] for x in swingingBallsol(t[(t .>= step3_start) .& (t .<= step3_end)])]

u = vcat(step1,step2,step3)
u = transpose(reduce(hcat,u))

using Javis
#using Luxor

step1 = transpose(reduce(hcat,step1))
θ1 = step1[:,1]
θ2 = step1[:,3]
θ3 = step1[:,5]

leg_outer_step1 = Point.(cos.(θ1)*L2 + cos.(θ2)*h0,sin.(θ1)*L2 + sin.(θ2)*h0)
leg_inner_step1 = Point.(cos.(θ1)*L2,sin.(θ1)*L2)
beam_tip_step1  = Point.(cos.(θ1)*0,sin.(θ1)*0)
beam_tail_step1 = Point.(cos.(θ1)*(L1+L2),sin.(θ1)*(L1+L2))
m2_outer_step1  = Point.(cos.(θ3)*L3,sin.(θ3)*L3)

step23 = transpose(reduce(hcat,vcat(step2,step3)))
θ1 = step23[:,1]
θ2 = step23[:,3]
θ3 = step23[:,5]

leg_outer_step23 = Point.(-cos.(θ2)*0  .+ rL                        ,-sin.(θ2)*0)
leg_inner_step23 = Point.(-cos.(θ2)*h0 .+ rL                        ,-sin.(θ2)*h0)
beam_tip_step23  = Point.(-cos.(θ2)*h0 .+ rL-cos.(θ1)*L2            ,-sin.(θ2)*h0-sin.(θ1)*L2)
beam_tail_step23 = Point.(-cos.(θ2)*h0 .+ rL+cos.(θ1)*L1            ,-sin.(θ2)*h0+sin.(θ1)*L1)
m2_outer_step23  = Point.(-cos.(θ2)*h0 .+ rL-cos.(θ1)*L2+cos.(θ3)*L3,-sin.(θ2)*h0-sin.(θ1)*L2+sin.(θ3)*L3)

leg_outer = vcat(leg_outer_step1,leg_outer_step23)
leg_inner = vcat(leg_inner_step1,leg_inner_step23)
beam_tip = vcat(beam_tip_step1,beam_tip_step23)
beam_tail = vcat(beam_tail_step1,beam_tail_step23)
m2_outer = vcat(m2_outer_step1,m2_outer_step23)

#print(leg_outer_step1)

function ground(args...)
    background("white") # canvas background
    sethue("black") # pen color
    scale(100,-100) # pixels per meter
end

function cap_line(p1,p2,color="black",thick=0.06)
    sethue(color)
    setline(thick*100)
    circle(p1,0.06, :fill)
    circle(p2,0.06, :fill)
    line(p1,p2, :stroke)
end

#scale(100) # pixels per meter

vid = Video(500,500)

Background(1:1000, ground)

#Luxor.scale(0.01,0.01) # meters per pixel
#Luxor.scale(100,100) # pixels per meter

beam = Object((video,object,frame) -> cap_line(beam_tip[frame],beam_tail[frame],"black"))
leg = Object((video,object,frame) -> cap_line(leg_outer[frame],leg_inner[frame],"gray"))
rope = Object((video,object,frame) -> cap_line(beam_tip[frame],m2_outer[frame],"orange",0.005))

#act!(beam,Action(1:1000,scale(100)))
#act!(leg,Action(1:1000,scale(100)))
#act!(rope,Action(1:1000,scale(100)))

render(
    vid,
    pathname="treb.gif",
    liveview=true
)

end
