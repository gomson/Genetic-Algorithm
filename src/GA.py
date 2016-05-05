# @Auther="Yuanhui Yang"
# @Email="yuanhui.yang@u.northwestern.edu"

from math import log, ceil, pow, exp, fabs
from random import random, randint, randrange

def fitnessfun(x): #输入popsize x UnitLength行向量
                global popsize
                Fitvalue = [0 for col in range(popsize)]
                for i in range(popsize):
                                temp = transform(x[i])
                                Fitvalue[i] = targetfun(temp)
                Fitvalue_Cumsum = Fitvalue[:]
                for i in range(1, popsize):
                                Fitvalue_Cumsum[i] =  Fitvalue_Cumsum[i-1] + Fitvalue[i]
                P_Cumsum = Fitvalue_Cumsum[:]
                for i in range(popsize):
                                P_Cumsum[i] /= Fitvalue_Cumsum[-1]
                y = [Fitvalue, P_Cumsum]
                return y#输出2 x popsize矩阵

def transform(x): #输入1 x UnitLength行向量, 即1条染色体
                global BitLength
                global boundsbegin
                global delta
                global VarNum
                y = [0 for col in range(VarNum)]
                for i in range(VarNum):
                                for j in range(BitLength[i]):
                                                y[i] += x[i * BitLength[i] + j] * pow(2, BitLength[i] - 1 - j)
                                y[i] *= delta[i]
                                y[i] += boundsbegin[i]
                return  y #输出1 x VarNum行向量, 即染色体对应自变量

def targetfun(x): #输入1 x VarNum行向量
                #边界处理
                global Initial
                g1 = 1.5 + x[0] * x[1] - x[0] - x[1]
                g2 = -x[0] * x[1]
                if g1 > 0 or g2 >10:
                                y = 0
                else:
                                y = Initial + exp(x[0]) * (4 * x[0] ** 2 + 2 * x[1] ** 2 + 4 * x[0] * x[1] + 2 * x[1] + 1) #\[\exp \left( {x\left[ 0 \right]} \right) \times \left( {4 \times x{{\left[ 0 \right]}^2} + 2 \times x{{\left[ 1 \right]}^2} + 4 \times x\left[ 0 \right] \times x\left[ 1 \right] + 2 \times x\left[ 1 \right] + 1} \right)\]
                                y = 1 / y
                return y #输出1个实数, 即染色体对应适应度值

def selection(x): #输入1 x UnitLength行向量
                y = [0 for col in range(2)]
                for i in range(2):
                                temp = random()
                                j = 0
                                while x[j] < temp:
                                                j += 1
                                y[i] = j
                return y #输出1 x 2行向量

def mutation(x): #输入1 x UnitLength行向量
                global pmutation
                global UnitLength
                y = x[:]
                temp1 = IfCroIfMut(pmutation)
                if temp1 == 1:
                                temp2 = randrange(UnitLength)
                                temp2 = int(temp2)
                                y[temp2] = fabs(x[temp2] - 1)
                return y #输出1 x UnitLength行向量

def IfCroIfMut(x): #输入1 x 1纯小数
                temp1 = [0 for col in range(10000)]
                temp2 = 10000 * x
                temp2 = ceil(temp2)
                temp2 = int(temp2)
                temp1[0: temp2] = [1 for col in range(temp2)]
                temp3 = randrange(10000)
                y = temp1[temp3]
                return y#输出1 x 1, 0或1

def crossover(x): #输入2 x UnitLength矩阵
                global UnitLength
                global pcrossover
                y = [[0 for col in range(UnitLength)] for row in range(2)]
                temp1 = IfCroIfMut(pcrossover)
                if temp1 == 1:
                                temp2 = randrange(UnitLength)
                                y[0] = x[0][0: temp2] + x[1][temp2: UnitLength]
                                y[1] = x[0][0: temp2] + x[0][temp2: UnitLength]
                else:
                                y[0] = x[0]
                                y[1] = x[1]
                return y #输出2 x UnitLength矩阵

def main():
                global BitLength #染色体长度(1维数组)
                global UnitLength #1条染色体长度
                global boundsbegin #变量下边界(矩形)
                global boundsend #变量上边界(矩形)
                global popsize #种群规模
                global pcrossover #交配概率
                global pmutation #变异概率
                global delta #编码精度
                global Initial #确保适应度为正值
                global VarNum #变量数
                Initial = 100 #确保适应度为正值
                pcrossover = 0.1 #交配概率
                pmutation = 0.4 #变异概率
                UnitLength = 0 #1条染色体长度
                boundsbegin = [-10, -10 ]#变量下边界(矩形)
                boundsend = [10, 10] #变量上边界(矩形)
                VarNum = len(boundsbegin) #变量数
                BitLength = [0 for col in range(VarNum)] #染色体长度(1维数组)
                delta = [0 for col in range(VarNum)] #编码精度
                precision = [1e-4, 1e-4]
                popsize = 100 #种群规模
                temp = popsize % 2
                if temp ==0:
                                popsize =popsize
                else:
                                popsize += 1
                Generationmax = 500
                for i in range(VarNum):
                                BitLength[i] = boundsend[i] - boundsbegin[i]
                                BitLength[i] /= precision[i]
                                BitLength[i] = log(BitLength[i])
                                BitLength[i] /= log(2)
                                BitLength[i] = ceil(BitLength[i])
                                BitLength[i] = int(BitLength[i])
                UnitLength = sum(BitLength)
                delta = BitLength[:]
                for i in range(VarNum):
                                delta[i] =  boundsend[i] - boundsbegin[i]
                                delta[i] /= pow(2, BitLength[i]) -1
                population = [[0 for col in range(UnitLength)] for row in range(popsize)]
                for i in range(popsize):
                                for j in range(UnitLength):
                                                temp = randint(0, 1)
                                                temp = int(temp)
                                                population[i][j] = temp
                Fitnessfun = fitnessfun(population)
                ValueMax = [0 for col in range(Generationmax)]
                NumMax = [0 for col in range(Generationmax)]
                VarMax = [0 for col in range(Generationmax)]
                temp4 = [[0 for col in range(UnitLength)] for row in range(popsize)]
                for i in range(Generationmax):
                                for j in range(0, popsize, 2):
                                                temp1 = selection(Fitnessfun[1])
                                                temp2 = [population[temp1[0]], population[temp1[1]]]
                                                temp3 = crossover(temp2)
                                                temp4[j] = mutation(temp3[0])
                                                temp4[j + 1] = mutation(temp3[1])
                                population = temp4[:]
                                Fitnessfun = fitnessfun(population)
                                temp = Fitnessfun[0]
                                ValueMax[i] = max(temp)
                                NumMax[i] = temp.index(ValueMax[i])
                                VarMax[i] = transform(population[NumMax[i]])
                print population
                print ValueMax
                print NumMax
                print VarMax
                print VarMax[-1]
main()
