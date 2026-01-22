import numpy as np 
import matplotlib.pyplot as plt 
from smooth import smooth_pde


def main():
    # set parameters
    sigma = 5. # gauss smoothing parameters, in km
    NGLL = 5 
    model_list = ['vp','vs','rho']
    model1d = 'prem'  # prem,ak135

    # smooth
    model,xx,xstore,displ,core_string = smooth_pde(NGLL,sigma,model_list,model1d)
    nglob = len(xstore)

    # write axisem external model
    f = open(f"{model1d}.smooth.bm","w")
    f.write(f"NAME {model1d}_smooth\n")
    f.write("ANELASTIC       T\n")
    f.write("ANISOTROPIC     T\n")
    f.write('UNITS           m\n')
    f.write("COLUMNS       radius      rho      vpv      vsv      qka      qmu      vph      vsh      eta\n")

    for iglob in range(nglob):
        vp0,vs0,rho0 = displ[iglob,:]
        f.write("\t %f %f %f %f %f %f %f %f %f\n"%(6371000 - 1000*xstore[iglob,0],rho0,vp0,vs0,
                                                9999.,9999.,vp0,vs0,1.))
    f.write(core_string)

    f = open(f"{model1d}.txt","w")
    for i in range(nglob):
        f.write("%f " %(-xstore[i,0]))
        for im in range(model.shape[-1]):
            f.write("%f " %(displ[i,im]))
        f.write("\n")
    f.close()

    plt.figure(1,figsize=(6,14))
    plt.plot(displ[:,1] / 1000,-xstore.flatten())
    plt.plot(model[:,:,1].flatten()/1000,-xx.flatten())
    plt.savefig("smooth.jpg")

if __name__ == "__main__":
    main()


