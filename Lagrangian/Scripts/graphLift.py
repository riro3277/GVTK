from matplotlib import pyplot as plt
import numpy as np

fig, axs = plt.subplots(2)
for p in range(0,9):
    path = "D:/lagrangiantest/PoiseuilleFlow/U1/ShearGradForceNoCorrection" + str(p) + ".txt"
    file = open(path, 'r')
    lines = file.readlines()
    sheargradmag = []
    sdfarr = []
    time = []
    for line in lines:
        temp = line.split(",")
        t = float(temp[0].split(":")[1])
        sdf = float(temp[1].split(":")[1])
        sheargrad = float(temp[2].split(":")[1])
        time.append(t)
        sdfarr.append(sdf)
        sheargradmag.append(sheargrad)

        print(len(time), len(sdfarr))
    axs[0].plot(time, sdfarr, label = "P="+str(p))
    axs[0].set_title("No Correction SDF vs Time")
    axs[0].set_ylim([0, 1])

    axs[1].plot(time, sheargradmag, label = "P="+str(p))
    axs[1].set_title("No Correction Shear Grad Magnitude vs Time")
    axs[1].set_ylim([0, 0.5])
axs[0].legend(loc='upper left')
axs[1].legend(loc='upper left')
plt.show()

    # plt.show()
    # # axs[0, p].plot(time, sdfarr)
    # axs[0, p].set_title('SDF vs Time p='+str(p))
    # axs[1, p].plot(time, sheargradmag, 'tab:orange')
    # axs[1, p].set_title('Shear Grad Magnitude vs Time p=' +str(p))
    #
    # for ax in axs.flat:
    #     ax.set(xlabel='x-label', ylabel='y-label')
    #
    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #     ax.label_outer()
