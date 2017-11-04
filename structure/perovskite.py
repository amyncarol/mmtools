import math

class StructureMap:
    #Shannon ionic radii with oxidation state and coordination  (and some with Zachariasen's relationship)
    A_radius = {'Li':1.13, 'Na':1.39, 'K':1.64, 'Rb':1.72, 'Cs':1.88, 'Tl':1.7}   #12 coodination

    B1_radius = {'Li':0.76, 'Na':1.02, 'K':1.38, 'Rb':1.52, 'Cs':1.67,              'Cu':0.77,'Ag':1.15,'Au':1.37, 'Tl':1.5, 'Hg': 1.19, 'In':1.4}  #6 coodination
           ####the 1+ radius for In was guessed based on In3+, Tl1+, Tl3+

    B2_radius = {'Bi':1.03,'Sb':0.76,'In':0.8,'Au':0.85,'Ta':0.72,'Mo':0.69,'Sc':0.745,'Y':0.9,'La':1.032,                 'Ce':1.01,'Pr':0.99,'Nd':0.983,'Sm':0.958,'Eu':0.947,'Gd':0.938,'Tb':0.923,'Dy':0.912,                 'Er':0.89,'Tm':0.88,'Lu':0.861, 'Cr':0.615, 'Fe':0.645, 'Ga':0.62, 'Ho':0.901, 'Yb':0.868,                 'Al':0.535, 'Tl':0.885, 'Rh':0.665, 'As':0.58, 'Ni':0.6, 'Co':0.61, 'Ir':0.68, 'Ru':0.68,                 'Mn':0.645, 'Nb':0.72, 'V':0.64, 'Ti':0.67, 'Pd':0.76 }   
            #6 coodination #Fe:high spin Ni:high spin Co: high spin Mn: high spin

    X_radius = {'F':1.33,'Cl':1.81,'Br':1.96,'I':2.2}   #2+4=6 coodination
    
    def get_tolerance_factor(self, A, B1, B2, X):
        try:
            return (self.A_radius[A]+self.X_radius[X])/((self.B1_radius[B1]+self.B2_radius[B2])/2+self.X_radius[X])/math.sqrt(2)
        except:
            return None
    
    def get_oct_factor(self, A, B1, B2, X):
        try:
            return (self.B1_radius[B1]+self.B2_radius[B2])/2/self.X_radius[X]
        except:
            return None

