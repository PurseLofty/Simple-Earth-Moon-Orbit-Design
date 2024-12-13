import math
import numpy as np

lamda_b = 50 #入口点经度
fai_b = 0 #入口点纬度
R_EA = 6578.14 #地球停泊轨道半径
i_L = 90 #绕月轨道倾角
epsilon = 4.243 #入口点倾角
R_EL = 384400 #地月距
V_L = 1.022 #月球公转线速度
miu_E = 398600 #地球引力参数
miu_L = 4902.8 #月球引力参数
rou = 66200 #月球影响球半径
angle = 23.44 #黄赤交角

lamda_b = math.radians(lamda_b)
fai_b = math.radians(fai_b)
i_L = math.radians(i_L)
epsilon = math.radians(epsilon) #全部转换为弧度

delta = math.acos(math.cos(lamda_b) * math.cos(fai_b))
pusai = math.asin(math.sin(fai_b) / math.sin(delta))
l = math.asin(math.sin(fai_b) / math.sin(i_L))

R_B = math.sqrt(R_EL**2 + rou**2 -2*R_EL*rou*math.cos(delta))
gama = math.asin(rou * math.sin(delta) / R_B)
d_0 = R_B * math.cos(gama)
e_0 = R_B * math.sin(gama) * math.cos(pusai)
f_0 = R_B * math.sin(gama) * math.sin(pusai)

def sgn(value):
    if value > 0:
        return 1
    elif value < 0:
        return -1
    else:
        return 0 #定义符号函数

if i_L < math.pi / 2:
    omiga_L = math.pi - lamda_b - math.acos(math.cos(l) / math.cos(fai_b))*sgn(fai_b)
else:
    omiga_L = math.pi - lamda_b + math.acos(math.cos(l) / math.cos(fai_b))*sgn(fai_b)

#入口点速度方向
matrix_RZ = np.array([[math.cos(-omiga_L), -math.sin(-omiga_L), 0],
                     [math.sin(-omiga_L), math.cos(-omiga_L), 0],
                     [0, 0, 1]])
matrix_iL = np.array([[1, 0, 0],
                     [0, math.cos(-i_L), -math.sin(-i_L)],
                     [0, math.sin(-i_L), math.cos(-i_L)]])
vector = np.array([[-math.cos(l - epsilon)],
                   [-math.sin(l - epsilon)],
                   [0]])
product_matrix = np.dot(matrix_RZ, matrix_iL)
v_B0Vector = np.dot(product_matrix, vector)
a_0, b_0, c_0 = v_B0Vector.flatten()

#入口点速度大小
m_0 = c_0 / a_0
n_0 = b_0 / a_0
g_0 = e_0 * V_L
h_0 = e_0*n_0 + d_0 + f_0*m_0
s_0 = 2 * V_L * n_0
p_0 = 1 + m_0**2 + n_0**2
t_0 = 2*g_0*h_0 - (R_B**2)*s_0
u_0 = g_0**2 - (R_B*V_L)**2
q_0 = h_0**2 - (R_B**2)*p_0
Q = 2 * miu_E * ((1/R_B) - (1/R_EA))
C_1 = (R_EA**2)*p_0 + q_0
C_2 = (R_EA**2)*s_0 + t_0
C_3 = u_0 - (R_EA**2)*(Q - V_L**2)
V_BX = (-C_2 - math.sqrt(C_2**2 - 4*C_1*C_3)) / (2*C_1)
V_BY = n_0*V_BX + V_L
V_BZ = m_0*V_BX
V_B = math.sqrt(V_BX**2 + V_BY**2 + V_BZ**2)

#地心轨道参数
V_EA = math.sqrt(V_B**2 - Q)
delta_1 = V_EA - math.sqrt(miu_E/R_EA)
niu_1 = R_EA*(V_EA**2) / miu_E
a_1 = R_EA / (2 - niu_1)
p_1 = R_EA * niu_1
e_1 = math.sqrt(1 - (p_1/a_1))
f_B = math.acos((1/e_1) * (p_1/R_B - 1))
E_B = 2*math.atan(math.sqrt((1-e_1)/(1+e_1)) * math.tan(f_B/2))

try:
    f_EA = math.acos((1/e_1) * (p_1/R_EA - 1))
except:
    f_EA = 0 #防止0附近的精度误差导致报错

E_EA = 2*math.atan(math.sqrt((1-e_1)/(1+e_1)) * math.tan(f_EA/2))
T_1 = math.sqrt(a_1**3 / miu_E) * ((E_B - e_1*math.sin(E_B)) - (E_EA - e_1*math.sin(E_EA)))
T_1 = T_1 / 3600

R_BVector = np.array([[d_0],
                      [e_0],
                      [f_0]])
scalar_1 = 1 - (R_EA*(1-math.cos(f_B-f_EA))/p_1)
V_BVector = np.array([[V_BX],
                      [V_BY],
                      [V_BZ]])
scalar_2 = R_EA * R_B * math.sin(f_B-f_EA) / math.sqrt(miu_E*p_1)
R_EAVector = (R_BVector * scalar_1) + (V_BVector * scalar_2)
x, y, z = R_EAVector.flatten()

#月心轨道参数
v_B = math.sqrt(p_0 * (V_BX**2))
niu_2 = rou * (v_B**2) / miu_L
p_2 = rou * niu_2 * (math.sin(epsilon)**2)
a_2 = rou / (2 - niu_2)
e_2 = math.sqrt(1 - (p_2/a_2))
f_B2 = math.acos((1/e_2) * (p_2/rou - 1))
H_B2 = 2 * math.atanh(math.sqrt((e_2-1)/(e_2+1)) * math.tan(f_B2/2))
T_2 = math.sqrt(-a_2**3 / miu_L) * (e_2*math.sinh(H_B2) - H_B2)
T_2 = T_2 / 3600
r_LP = p_2 / (1 + e_2)
v_LP = math.sqrt((v_B**2) - 2*miu_L*((1/rou) - (1/r_LP)))
delta_2 = v_LP - math.sqrt(miu_L/r_LP)

T = T_1 + T_2
energy_cost = delta_1 + delta_2

print('在地球停泊轨道上变轨处的地心位置矢量为:', x, y, z)
print('总时间为（小时）：', T)
print('近月距为（千米）：', r_LP)
print('两次变轨所需能量为：', energy_cost)