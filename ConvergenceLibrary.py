"""
Library of function declearations used in plotting verificaiton results
"""
import scipy as sp
import matplotlib
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import cm
from scipy import linalg
from shutil import copyfile
import os 
import sys



 
class ThermoData(object):
    Temperature=False;
    Density=False;
    Enthalpy=False;
    StaticEnergy=False;
    Entropy=False;
    Pressure=False;
    SoundSpeed2=False;
    dPdrho_e=False;
    dPde_rho=False;
    dTdrho_e=False;
    dTde_rho=False;
    Cp=False;
    Mu=False;
    Kt=False;
    
    def load_data(self, filename):
        self.data = sp.genfromtxt(filename, delimiter=', ')
        self.Temperature  = self.data[:,0];
        self.Density      = self.data[:,1];
        self.Enthalpy     = self.data[:,2];
        self.StaticEnergy = self.data[:,3];
        self.Entropy      = self.data[:,4];
        self.Pressure     = self.data[:,5];
        self.SoundSpeed2  = self.data[:,6];
        self.dPdrho_e     = self.data[:,7];
        self.dPde_rho     = self.data[:,8];
        self.dTdrho_e     = self.data[:,9];
        self.dTde_rho     = self.data[:,10];
        self.Cp           = self.data[:,11];
        self.Mu           = self.data[:,12];
        #self.dmudrho_T    = self.data[:,13];
        #self.dmudT_rho    = self.data[:,14];
        self.Kt           = self.data[:,15];
        #self.dktdrho_T    = self.data[:,16];
        #self.dktdT_rho    = self.data[:,17];
        return 
        
    def load_FP_output(self, filename):
        randoms = sp.genfromtxt(filename, skip_header=2)
        self.Density      = randoms[:,0];
        self.Pressure     = randoms[:,1];
        self.SoundSpeed2  = randoms[:,2]**2;
        self.Cp           = randoms[:,3];
        self.Entropy      = randoms[:,4];
        self.Mu           = randoms[:,5];
        self.Kt           = randoms[:,6];
        self.dPdrho_e     = randoms[:,7];
        self.dPde_rho     = randoms[:,8];
        self.dTdrho_e     = randoms[:,9];
        self.dTde_rho     = randoms[:,10];
        self.Temperature  = randoms[:,11];
        self.StaticEnergy = randoms[:,12];
        self.Enthalpy     = randoms[:,13];
        return
        

class RandomSamples(ThermoData):
    
    def __init__(self,filename):
        print 'Loading random verification data'
        self.load_FP_output(filename);
        print 'DONE Loading random verification data'
        
        #Prepare the input files
        print 'Preparing input files...'
        sp.savetxt('rhoe_in.dat',sp.column_stack((self.Density,self.StaticEnergy)), delimiter='\t')
        sp.savetxt('PT_in.dat'  ,sp.column_stack((self.Pressure,self.Temperature)), delimiter='\t')
        sp.savetxt("Prho_in.dat",sp.column_stack((self.Pressure,self.Density)), delimiter='\t')
        sp.savetxt("rhoT_in.dat",sp.column_stack((self.Density,self.Temperature)), delimiter='\t')
        sp.savetxt("Ps_in.dat"  ,sp.column_stack((self.Pressure,self.Entropy)), delimiter='\t')
        sp.savetxt("hs_in.dat"  ,sp.column_stack((self.Enthalpy,self.Entropy)), delimiter='\t')        
        print 'DONE Preparing input files...'
        return


class RestartSamples(object):
    Density=None;
    Pressure=None;
    Temperature=None;
    restarts = None;
    
    def __init__(self,filename, dim='2D'):
        restarts = sp.genfromtxt(filename, skip_header=2)
        if dim == '2D':
            self.Density      = restarts[:,3];
            self.Pressure     = restarts[:,9];
            self.Temperature  = restarts[:,10];        
        elif dim == '3D':
            self.Density      = restarts[:,4];
            self.Pressure     = restarts[:,9];
            self.Temperature  = restarts[:,10];        
        self.restarts=restarts;
            
        print 'Restart points: ', len(self.Density)
        return;
        
    def plot_restarts(self, scatter_x ='Density', scatter_y='Pressure'):
        plt.scatter(getattr(self,scatter_x), getattr(self,scatter_y));
        return;


        
class SciPy_InterpolatedData(ThermoData):
    

    def __init__(self, which_case, LUT, RandomSamples, interp_type):
        print 'SciPy Interpolating ', which_case
        
        select = {\
        "rhoe":('Density','StaticEnergy'),\
        "PT":('Pressure','Temperature'),\
        "Prho":('Pressure','Density'),\
        "rhoT":('Density','Temperature'),\
        "Ps":('Pressure','Entropy'),\
        "hs":('Enthalpy','Entropy')\
        }
        
        thermo1, thermo2, = select[which_case]
        x =getattr(LUT,thermo1)
        y =getattr(LUT,thermo2)
        samples_x = getattr(RandomSamples,thermo1)
        samples_y = getattr(RandomSamples,thermo2)
        setattr(self,thermo1, samples_x)
        setattr(self,thermo2, samples_y)
        
        variables = sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
        'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho',\
        'dTdrho_e','dTde_rho','Cp','Mu','Kt']);
        
        for var in variables[sp.where((variables!=thermo1) * (variables!=thermo2))]:
            z = getattr(LUT,var)            
            interp_func = sp.interpolate.griddata((x,y),z,sp.column_stack((samples_x,samples_y)),\
            method=interp_type) 
            nan_index = sp.where(sp.isnan(interp_func))
            interp_func[nan_index]= sp.interpolate.griddata((x,y),z,\
            sp.column_stack((samples_x[nan_index],samples_y[nan_index])),\
            method='nearest') 
            setattr(self,var,interp_func)
            
        return  
          
        
class SU2_InterpolatedData(ThermoData):
    name='';
    data=False;
    
    def __init__(self, filename):
        self.load_data(filename);

        
class PRGrid(ThermoData):
    
    meshfile=False;
    P_dim=0;
    D_dim=0;
        
    def __init__(self, filename):
        self.load_data(filename);
        self.P_dim = sp.where((self.Density-self.Density[0])!=0)[0][0]
        self.D_dim = len(self.Density)/self.P_dim
        print 'P-RHO dimensions: %i by %i'%(self.P_dim, self.D_dim);

    def plot_mesh(self,thermo_x,thermo_y):
        x = getattr(self,thermo_x)
        y = getattr(self,thermo_y)
        
        for i in range(self.P_dim):
            plt.plot(x[i:][::self.P_dim], y[i:][::self.P_dim], 'c-', alpha=0.3);
        for i in range(self.D_dim):
            plt.plot(x[i*self.P_dim:(i+1)*self.P_dim], y[i*self.P_dim:(i+1)*self.P_dim], 'c-', alpha=0.3);

class PRGridSkewed(ThermoData):
    
    meshfile=False;
    P_dim=0;
    D_dim=0;
        
    def __init__(self, filename):
        self.load_FP_output(filename)
        self.P_dim = sp.where((self.Density-self.Density[0])!=0)[0][0]
        self.D_dim = len(self.Density)/self.P_dim
        print 'P-RHO dimensions: %i by %i'%(self.P_dim, self.D_dim);

    def plot_mesh(self,thermo_x,thermo_y):
        x = getattr(self,thermo_x)
        y = getattr(self,thermo_y)
        
        for i in range(self.P_dim):
            plt.plot(x[i:][::self.P_dim], y[i:][::self.P_dim], 'c-', alpha=0.3);
        for i in range(self.D_dim):
            plt.plot(x[i*self.P_dim:(i+1)*self.P_dim], y[i*self.P_dim:(i+1)*self.P_dim], 'c-', alpha=0.3);

class RefinementLevel(object):
   cases=False;
   filename=False;
   variables=sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
        'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho',\
        'dTdrho_e','dTde_rho','Cp','Mu','Kt']);
   select = {\
        "rhoe":('Density','StaticEnergy'),\
        "PT":('Pressure','Temperature'),\
        "Prho":('Pressure','Density'),\
        "rhoT":('Density','Temperature'),\
        "Ps":('Pressure','Entropy'),\
        "hs":('Enthalpy','Entropy')\
        }
   LUT={};
   SU2={};
   RandomSamples=False;
   SciPy={}; 
            
   def __init__(self, filename, cases):        
       self.cases =cases;        
       self.filename=filename;
     
   def load_random_samples(self,RandomSamples):
       self.RandomSamples = RandomSamples
   
   def load_mesh(self):
       self.LUT = PRGrid('mesh.dat')
       
   def load_results_SU2(self):
       self.SU2={i: SU2_InterpolatedData(str(i+'_out.dat')) for i in self.cases};            
            
   def load_results_SciPy(self,interp_kind='linear'):
       self.SciPy={i:SciPy_InterpolatedData(i,self.LUT, self.RandomSamples, \
       interp_kind ) for i in self.cases};
   
   def plot_hist_compare(self,which_case):
        plt.ylabel('Percentage of points')
        plt.xlabel('Percentage RMS relative error')
        
        def yto_percent(y, x):
            s = str(sp.around((y/(len(self.REL_ERR)*1.0)*100),2))
            if matplotlib.rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                return s + '%'     

        def xto_percent(y, x):
            s = str(y*100)
            if matplotlib.rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                    return s + '%' 
        
        thermo1, thermo2, = self.select[which_case]
        #Plot the SU2 error
        i=0;
        self.REL_ERR = 0;
        for v in self.variables[sp.where\
        ((self.variables!=thermo1) * (self.variables!=thermo2))]:
            i=i+1;
            self.REL_ERR = self.REL_ERR + \
            ((getattr(self.SU2[which_case],v)-getattr(self.RandomSamples,v))/\
            (getattr(self.RandomSamples,v)))**2;
        self.REL_ERR = sp.sqrt(self.REL_ERR/i)
        plt.hist(self.REL_ERR, bins=25, color='k', alpha=0.3, label='SU2')
        print 'Error max SU2', max(self.REL_ERR)
        
        #Plot the SciPy error
        i =0;
        self.REL_ERR = 0;
        for v in self.variables[sp.where\
        ((self.variables!=thermo1) * (self.variables!=thermo2))]:
            i=i+1;
            self.REL_ERR = self.REL_ERR + \
            ((getattr(self.SciPy[which_case],v)-getattr(self.RandomSamples,v))/\
            (getattr(self.RandomSamples,v)))**2;
        self.REL_ERR = sp.sqrt(self.REL_ERR/i)
        
        plt.hist(self.REL_ERR, bins=25, color='c', alpha=0.5, label='SciPy')
        print 'Error max SciPy', max(self.REL_ERR)

        
        formatter_y = FuncFormatter(yto_percent)
        formatter_x = FuncFormatter(xto_percent)
        plt.gca().yaxis.set_major_formatter(formatter_y)
        plt.gca().xaxis.set_major_formatter(formatter_x)
        plt.grid(which='both')
        plt.legend()

       
        return
       
              
if __name__ == '__main__':
    PR_files = os.listdir('TableLibrary')
    levels = len(PR_files);
    cases = ["rhoe","PT", "Prho","rhoT","Ps","hs"]
    Refinement_Levels = [RefinementLevel(PR_files[i], cases) for i in range(levels)]
    Random_Samples = RandomSamples('random.dat')
    for level in Refinement_Levels:
        output_files = ["rhoe_out.dat","PT_out.dat", "Prho_out.dat","rhoT_out.dat","Ps_out.dat","hs_out.dat"]
        for f in output_files:
            try:
                 os.remove(f)
            except OSError:
                 pass
             
        #Use the thermotable related to the desired refinement level
        copyfile('TableLibrary/'+level.filename, 'CO2.rgp')  
        print 'Mesh file:',  level.filename;
        #Recompile if necessary             
        os.system('make Debug/makefile');
        print 'Running...'
        #Run and save the log
        os.system('Debug/TableReader>log');
        print 'DONE'
        level.load_mesh();
        level.load_random_samples(Random_Samples)
        level.load_results_SU2();
        level.load_results_SciPy('linear');
        #level.LUT.plot_mesh('Density','StaticEnergy')
        plt.figure(figsize=(16,6))
        level.plot_hist_compare('rhoe')
        #level.plot_hist_compare('hs')
        #level.plot_hist_compare('Prho')
"""           


def plot_hist_RMS(REL_ERR):
    #adapted from http://matplotlib.org/examples/pylab_examples/histogram_percent_demo.html
    plt.figure(figsize=(16,5))
    plt.title('Historgram of RMS Relative Error')
    plt.ylabel('Percentage of points')
    plt.xlabel('Percentage RMS relative error')
    plt.hist(REL_ERR, bins=25, color='g', alpha=0.5)
    formatter_y = FuncFormatter(yto_percent)
    formatter_x = FuncFormatter(xto_percent)
    plt.gca().yaxis.set_major_formatter(formatter_y)
    plt.gca().xaxis.set_major_formatter(formatter_x)
    plt.grid(which='both')
    return

def yto_percent(y, x):
        s = str(sp.around((y/(len(REL_ERR)*1.0)*100),2))
        if matplotlib.rcParams['text.usetex'] is True:
            return s + r'$\%$'
        else:
            return s + '%' 
    

def xto_percent(y, x):
        s = str(y*100)
        if matplotlib.rcParams['text.usetex'] is True:
            return s + r'$\%$'
        else:
            return s + '%' 
        

def plot_hist_Individual(input1, input2):
    names = [["Temperature", interp_Temperature,T],\
            ["Density", interp_Density,rho],\
            ["Enthalpy",interp_Enthalpy,Enthalpy],\
            ["StaticEnergy", interp_StaticEnergy,StaticEnergy],\
            ["Entropy",interp_Entropy,Entropy],\
            ["Pressure",interp_Pressure,P],\
            ["SoundSpeed", interp_SoundSpeed2,SoundSpeed],\
            ["dPdrho_e",interp_dPdrho_e,dPdrho_e],\
            ["dPde_rho", interp_dPde_rho,dPde_rho],\
            ["Tdrho_e",interp_dTdrho_e,dTdrho_e],\
            ["dTde_rho",interp_dTde_rho,dTde_rho],\
            ["Cp",interp_Cp,Cp],\
            ["Mu", interp_Mu,Mu],\
            ["Kt", interp_Kt,Kt]]
    i =1
    plt.figure(figsize=(16,20))
    for n in names:
        rel_err = (n[1] - n[2])/max(n[2])
        plt.subplot(7,2,i)
        if (n[0] == input1) or (n[0] == input2):
            col ='b'
            plt.title(n[0]+" (input)")
        else:
            col ='g'
            plt.title(n[0])
        plt.hist(rel_err, bins=25, normed=False, color=col, alpha=0.5)
        formatter_y = FuncFormatter(yto_percent)
        formatter_x = FuncFormatter(xto_percent)
        plt.gca().yaxis.set_major_formatter(formatter_y)
        plt.gca().xaxis.set_major_formatter(formatter_x)
        plt.grid(which='both')
        
        i = i+1

    return




#Plot the interpolation grid: 

        
"""           

