import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.integrate import ode
class AD(object):
    def __init__(self, left, right, product, n):
        self.l = left
        self.r = right
        self.p = product
        self.n = n
        self.Rmat = np.array([['-' + self.n + 'kf', self.n + 'kb'],
                              ['-' + self.n + 'kf', self.n + 'kb'],
                              [self.n + 'kf', '-' + self.n + 'kb']],dtype=object)
        self.Rcnt = np.array([[self.l + '*' + self.r], [self.p]],dtype=object)
        self.Rxn = np.array([[self.l], [self.r], [self.p]],dtype=object)
class ES(object):
     def __init__(self, enzyme, substrate, product,n):
         self.enz = enzyme
         self.sub = substrate
         self.product = product
         self.n = n
         self.Rmat = np.array([['-'+self.n+'kf', self.n+'kb'],
                           ['-'+self.n+'kf', self.n+'kb' + '+' +self.n + 'kr'],
                           [self.n+'kf', '-'+self.n+'kb' +'+' +'-' +self.n + 'kr'],
                           ['', self.n+'kr']],dtype=object)

         self.Rcnt = np.array([[self.enz + '*'+self.sub],[self.enz+'-'+self.sub]],dtype=object)
         self.Rxn  = np.array([[self.enz],[self.sub],[self.enz+'-'+self.sub],[self.product]],dtype=object)
class model(object):
    def __init__(self,a,b):
        self.a = a
        self.b = b
        self.Rxn  = np.vstack((self.a.Rxn,self.b.Rxn))
        self.Rcnt = np.vstack((self.a.Rcnt,self.b.Rcnt))
        c = np.array(np.tile('', (self.a.Rmat.shape[0] + self.b.Rmat.shape[0], self.a.Rmat.shape[1] + self.b.Rmat.shape[1])),
                     dtype=object)
        n = self.a.Rmat.shape[0]
        m = self.a.Rmat.shape[1]
        c[0:n, 0:m] = self.a.Rmat
        c[n:, m:] = self.b.Rmat
        self.Rmat = c
    def fuse(self):
        delete = []
        Delete = []
        n=self.a.Rxn.shape[0]
        m=self.Rxn.shape[0]
        for i in range(n,m):
            for j in range(0,self.a.Rxn.shape[0]):
                if (self.Rxn[i] ==self.Rxn[j]):
                    self.Rmat[j,0:]=self.Rmat[j,0:]+self.Rmat[i,0:]
                    delete=delete+[i]
        self.Rxn = np.delete(self.Rxn,delete,axis=0)
        self.Rmat= np.delete(self.Rmat,delete,axis=0)
        for i in range(self.a.Rcnt.shape[0],self.Rcnt.shape[0]):
            for j in range(0,self.a.Rcnt.shape[0]):
                if (self.Rcnt[i] ==self.Rcnt[j]):
                    self.Rmat[0:,j]=self.Rmat[0:,j]+self.Rmat[0:,i]
                    Delete = Delete + i
        self.Rcnt=np.delete(self.Rcnt,Delete,axis=0)
        self.Rmat=np.delete(self.Rmat,Delete,axis=1)
source = AD('A', 'P', 'AP', 'R1')                #reactant_1,reactant_2, Product
gate = AD('P', 'E', 'PE', 'R2')                  #reactant_1,reactant_2, Product
drain = ES('A', 'E', '~A', 'R3')                 #enzyme,substrate,product
c = model(source, gate)
c.fuse()
c = model(c,drain)
c.fuse()
#print (c.Rcnt)
print (c.Rxn)
s=10 #A conc.
#  #
#print (c.Rmat)
i_p = 23.33
i_a = 45.95
i_e = 9.69
init_A = 0.33
init_P = 0.33
init_E = 49.45
Names_Rxn = c.Rxn
Names_Rcnt = c.Rcnt
Rmat_dict = {}
#all parameters initiated to zero
for i in range(0, c.Rmat.shape[0]):
        for j in range(0, c.Rmat.shape[1]):
            temp = c.Rmat[i, j]
            if ('+' in temp):
                temp1 = temp.split('+')[0]
                temp2 = temp.split('+')[1]
                if ('-' in temp1):
                    temp1 = temp1.split('-')[1]
                elif ('-' in temp2):
                    temp2 = temp2.split('-')[1]
                if ~(temp1 in Rmat_dict.keys()):
                    Rmat_dict[temp1] = 0
                elif ~(temp2 in Rmat_dict.keys()):
                    Rmat_dict[temp2] = 0
            elif ('-' in temp):
                if ~(temp.split('-')[1] in Rmat_dict.keys()):
                    Rmat_dict[temp.split('-')[1]] = 0
            elif ~(temp in Rmat_dict.keys()):
                Rmat_dict[temp] = 0
    # Original Matrix formula -> c.Rxn = c.Rmat*c.Rcnt
    # dictionary form         -> Rxn_dict,Rmat_dict, Rcnt_dict
#print (Rmat_dict.keys())
def ds_dt_ramp(t,S,P):
    global i_a
    global i_p
    global i_e
    global Rmat_dict
    Rmat_dict[""] = P[0]
    Rmat_dict["R1kf"] = P[1]
    Rmat_dict["R1kb"] = P[2]
    Rmat_dict["R2kf"] = P[3]
    Rmat_dict["R2kb"] = P[4]
    Rmat_dict["R3kf"] = P[5]
    Rmat_dict["R3kb"] = P[6]
    Rmat_dict["R3kr"] = P[7]
    global c
    global Names_Rcnt
    global Names_Rxn
    MAT = np.zeros_like(c.Rmat)
    # formula parsing MAT
    for i in range(0, c.Rmat.shape[0]):
        for j in range(0, c.Rmat.shape[1]):
            temp = c.Rmat[i, j]
            if ('+' in temp):
                temp1 = temp.split('+')[0]
                temp2 = temp.split('+')[1]
                if ('-' in temp1):
                    val1 = -1 * Rmat_dict[temp1.split('-')[1]]
                else:
                    val1 = Rmat_dict[temp1]
                if ('-' in temp2):
                    val2 = -1 * Rmat_dict[temp2.split('-')[1]]
                else:
                    val2 = Rmat_dict[temp2]
                MAT[i, j] = val1 + val2
            elif ('-' in temp):
                MAT[i, j] = -1 * Rmat_dict[temp.split('-')[1]]
            else:
                MAT[i, j] = Rmat_dict[temp]
    # Rmat_dict has an unordered parameter mapping
    # params is the ordered set of parameter names
    # params_value is the ordered set of parameter values
    # formula parsing RXN
    rxn_values = np.reshape(S,(7,))  #### This is where the S_vec goes,There's an adhoc Jugad of size changing here
    Rxn_dict = dict(zip([i[0] for i in Names_Rxn], rxn_values))
    rcnt_values = np.zeros_like(Names_Rcnt)
    # formula parsing RCNT
    for i in range(0, Names_Rcnt.shape[0]):
        if ('*' in Names_Rcnt[i, 0]):
            rcnt_values[i, 0] = Rxn_dict[Names_Rcnt[i, 0].split('*')[0]] * Rxn_dict[Names_Rcnt[i, 0].split('*')[1]]
        else:
            rcnt_values[i, 0] = Rxn_dict[Names_Rcnt[i, 0]]
    #print(rxn_values)
    #print(rcnt_values)
    #print(MAT)
    Rcnt_dict = dict(zip([i[0] for i in Names_Rcnt], rcnt_values))
    ### dS_dt is being defined
    dS_dt = np.dot(MAT,rcnt_values) + np.array([[i_a],[20*t],[0],[i_e],[0],[0],[0]])
    return dS_dt
r = ode(ds_dt_ramp)
r.set_integrator("vode")
r.set_initial_value(np.array([[init_A],[init_P],[0],[init_E],[0],[0],[0]]))
r.set_f_params(np.array([0.0,9.5,0.05,48.67,9.74,9.89,0.02,95.00]))
solution_ramp =np.vstack((np.array([0]),np.array([[init_A],[init_P],[0],[init_E],[0],[0],[0]])))
#print (solution)
T_max = 10.0
dt = 0.01
#print (r.integrate(r.t+dt))
while r.successful() and r.t < T_max-dt:
     a = np.array([r.t+dt])
     b = r.integrate(r.t+dt)
     solution_ramp=np.hstack((solution_ramp,np.vstack((a,b))))
#error = 2*np.linspace(0.0,3.0,31)-solution[1,:]
#print (parameter_names_in_order)
#print(solution[0,:])
#print(solution[0,:])
#print (np.sum(error**2))

# def ds_dt(t,S,P):
#     global i_a
#     global i_p
#     global i_e
#     global Rmat_dict
#     Rmat_dict[""] = P[0]
#     Rmat_dict["R1kf"] = P[1]
#     Rmat_dict["R1kb"] = P[2]
#     Rmat_dict["R2kf"] = P[3]
#     Rmat_dict["R2kb"] = P[4]
#     Rmat_dict["R3kf"] = P[5]
#     Rmat_dict["R3kb"] = P[6]
#     Rmat_dict["R3kr"] = P[7]
#     global c
#     global Names_Rcnt
#     global Names_Rxn
#     MAT = np.zeros_like(c.Rmat)
#     # formula parsing MAT
#     for i in range(0, c.Rmat.shape[0]):
#         for j in range(0, c.Rmat.shape[1]):
#             temp = c.Rmat[i, j]
#             if ('+' in temp):
#                 temp1 = temp.split('+')[0]
#                 temp2 = temp.split('+')[1]
#                 if ('-' in temp1):
#                     val1 = -1 * Rmat_dict[temp1.split('-')[1]]
#                 else:
#                     val1 = Rmat_dict[temp1]
#                 if ('-' in temp2):
#                     val2 = -1 * Rmat_dict[temp2.split('-')[1]]
#                 else:
#                     val2 = Rmat_dict[temp2]
#                 MAT[i, j] = val1 + val2
#             elif ('-' in temp):
#                 MAT[i, j] = -1 * Rmat_dict[temp.split('-')[1]]
#             else:
#                 MAT[i, j] = Rmat_dict[temp]
#     # Rmat_dict has an unordered parameter mapping
#     # params is the ordered set of parameter names
#     # params_value is the ordered set of parameter values
#     # formula parsing RXN
#     rxn_values = np.reshape(S,(7,))  #### This is where the S_vec goes,There's an adhoc Jugad of size changing here
#     Rxn_dict = dict(zip([i[0] for i in Names_Rxn], rxn_values))
#     rcnt_values = np.zeros_like(Names_Rcnt)
#     # formula parsing RCNT
#     for i in range(0, Names_Rcnt.shape[0]):
#         if ('*' in Names_Rcnt[i, 0]):
#             rcnt_values[i, 0] = Rxn_dict[Names_Rcnt[i, 0].split('*')[0]] * Rxn_dict[Names_Rcnt[i, 0].split('*')[1]]
#         else:
#             rcnt_values[i, 0] = Rxn_dict[Names_Rcnt[i, 0]]
#     #print(rxn_values)
#     #print(rcnt_values)
#     #print(MAT)
#     Rcnt_dict = dict(zip([i[0] for i in Names_Rcnt], rcnt_values))
#     ### dS_dt is being defined
#     dS_dt = np.dot(MAT,rcnt_values) + np.array([[i_a],[i_p],[0],[i_e],[0],[0],[0]])
#     return dS_dt
# r = ode(ds_dt)
# r.set_integrator("vode")
# r.set_initial_value(np.array([[init_A],[init_P],[0],[init_E],[0],[0],[0]]))
# r.set_f_params(np.array([0.0,9.5,0.05,9.5,0.01,9.89,0.02,95.00]))
# solution = np.vstack((np.array([0]),np.array([[init_A],[init_P],[0],[init_E],[0],[0],[0]])))
# #print (solution)
# T_max = 10.0
# dt = 0.01
# #print (r.integrate(r.t+dt))
# while r.successful() and r.t < T_max-dt:
#      a = np.array([r.t+dt])
#      b = r.integrate(r.t+dt)
#      solution=np.hstack((solution,np.vstack((a,b))))
# #error = 2*np.linspace(0.0,3.0,31)-solution[1,:]
# #print (parameter_names_in_order)
# #print(solution[0,:])
# #print(solution[0,:])
# #print (np.sum(error**2))
# currents = np.vstack((solution[0, 0:999], 100 * (solution[1:, 1:1000] - solution[1:, 0:999]))) #forward difference derivative 300 data points
currents_ramp = np.vstack((solution_ramp[0,0:999],100*(solution_ramp[1:,1:1000]-solution_ramp[1:,0:999])))
fig = plt.figure()
fig.suptitle('Ramp Response', fontweight='bold',fontsize=12)
graph = fig.add_subplot(111)
graph.set_title('IO Charecteristic')
graph.set_xlabel(r'Input Current $i_p[\frac{microMoles}{litre.sec}]$')
graph.set_ylabel(r'Output Current $i_{\~A}[\frac{microMoles}{litre.sec}]$')
fig.text(0.75,0.03, r'[$\frac{microMoles}{(litre.sec)}$]',fontsize=14)
fig.text(0.2,0.65,r'$y=x$',fontsize=15)
fig.text(0.02,0.3, 'Rate Constants \n [Approriate Dimensions]',fontsize=12,fontstyle='italic')
plt.subplots_adjust(left=0.1, top=0.9, bottom=0.4)
#l, = plt.plot(solution[0,:],solution[7,:]) ##Here t varies as 0,0.01, 0.02, ... 3.0 301 in all
#l, = plt.plot(currents[0,:],currents[7,:])
l,= plt.plot(20*currents_ramp[0,:],currents_ramp[7,:])
plt.plot(5*currents_ramp[0,:],50-5*currents_ramp[0,:])
plt.axis([0,100,0,50])       # RRR range to be shown on the axes x in (0,1) y in (-10,10)
axcolor = 'lightgoldenrodyellow'
axkf1 = plt.axes([0.1, 0.25, 0.2, 0.03], axisbg=axcolor) #in fraction of the screen-[left horizontal margin, lower vertical margin,length of the bar,thickness]
axkb1 = plt.axes([0.1, 0.21, 0.2, 0.03], axisbg=axcolor)
axkf2 = plt.axes([0.1, 0.17, 0.2, 0.03], axisbg=axcolor)
axkb2 = plt.axes([0.1, 0.13, 0.2, 0.03], axisbg=axcolor)
axkf3 = plt.axes([0.1, 0.09, 0.2, 0.03], axisbg=axcolor)
axkb3 = plt.axes([0.1, 0.05, 0.2, 0.03], axisbg=axcolor)
axkr3 = plt.axes([0.1, 0.01, 0.2, 0.03], axisbg=axcolor)

axi_a = plt.axes([0.65, 0.25, 0.2, 0.03], axisbg=axcolor) # 0.06,0.95,0.2,0.03
#axi_p = plt.axes([0.06, 0.91, 0.2, 0.03], axisbg=axcolor)
axi_e = plt.axes([0.65, 0.21, 0.2, 0.03], axisbg=axcolor) # 0.06,0.91,0.2,0.03

axi_init_A = plt.axes([0.65, 0.17, 0.2, 0.03], axisbg=axcolor)
axi_init_P = plt.axes([0.65, 0.13, 0.2, 0.03], axisbg=axcolor)
axi_init_E = plt.axes([0.65, 0.09, 0.2, 0.03], axisbg=axcolor)

skf1 = Slider(axkf1, 'Ks', 0.1, 10.0, valinit=9.5)
skb1 = Slider(axkb1, '`Ks', 0.001, 0.1, valinit=0.05)
skf2 = Slider(axkf2, 'Kg', 0.1, 50.0, valinit=48.67)
skb2 = Slider(axkb2, '`Kg', 0.1, 10.0, valinit=9.74)
skf3 = Slider(axkf3, 'Kd', 0.1, 10.0, valinit=9.89)
skb3 = Slider(axkb3, '`Kd', 0.001, 0.1, valinit=0.02)
skr3 = Slider(axkr3, 'Kr', 0.1, 100.0, valinit=95)
si_a = Slider(axi_a, 'i_a', 0.1, 50.0, valinit=45.95)
#si_p = Slider(axi_p, 'i_p', 0.1, 50.0, valinit=23.33)
si_e = Slider(axi_e, 'i_e',0.1,50.0,valinit=9.69)
s_init_A = Slider(axi_init_A, 'init_A',0.1,50.0,valinit=0.33)
s_init_P = Slider(axi_init_P, 'init_P',0.1,50.0,valinit=0.33)
s_init_E = Slider(axi_init_E, 'init_E',0.1,100.0,valinit=49.45)

def update(val):
    kf1 = skf1.val
    kb1 = skb1.val
    kf2 = skf2.val
    kb2 = skb2.val
    kf3 = skf3.val
    kb3 = skb3.val
    kr3 = skr3.val
    global ds_dt
    global i_a
    global i_p
    global i_e
    global init_A
    global init_P
    global init_E
    i_a = si_a.val
    i_p = si_p.val
    i_e = si_e.val
    init_A = s_init_A.val
    init_P = s_init_P.val
    init_E = s_init_E.val
    r = ode(ds_dt_ramp)
    r.set_integrator("vode")
    r.set_initial_value(np.array([[init_A], [init_P], [0], [init_E], [0], [0], [0]]))
    r.set_f_params(np.array([0.0, kf1,kb1,kf2,kb2,kf3,kb3,kr3]))
    solution = np.vstack((np.array([0]), np.array([[init_A], [init_P], [0], [init_E], [0], [0], [0]])))
    # print (solution)
    T_max = 10.0
    dt = 0.01
    # print (r.integrate(r.t+dt))
    while r.successful() and r.t < T_max - dt:
        a = np.array([r.t + dt])
        b = r.integrate(r.t + dt)
        solution = np.hstack((solution, np.vstack((a, b))))
    # error = 2*np.linspace(0.0,3.0,31)-solution[1,:]
    # print (parameter_names_in_order)
    # print(solution[0,:])
    # print(solution[0,:])
    # print (np.sum(error**2))
    currents_ramp = np.vstack((solution[0,0:999],100*(solution[1:, 1:1000] - solution[1:, 0:999])))  # forward difference derivative 300 data points
#    l.set_ydata(solution[7,:])      ##Here
    l.set_ydata(currents_ramp[7,:])       ##Here
    fig.canvas.draw_idle()
skf1.on_changed(update)
skb1.on_changed(update)
skf2.on_changed(update)
skb2.on_changed(update)
skf3.on_changed(update)
skb3.on_changed(update)
skr3.on_changed(update)
si_a.on_changed(update)
#si_p.on_changed(update)
si_e.on_changed(update)
s_init_A.on_changed(update)
s_init_E.on_changed(update)
s_init_P.on_changed(update)
# def reset(event):
#     sfreq.reset()
#     samp.reset()
plt.show()