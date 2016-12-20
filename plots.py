import math
import numpy as np
import matplotlib.pyplot as plt
import datetime


#plot time against height
def plot_th(dates,final_h):
    plt.plot(dates, final_h) 
    plt.gcf().autofmt_xdate()
    #range_x = max(final_time)-min(final_time)
    #range_y = max(final_h)-min(final_h)
    #axes.set_xlim([min(final_time)-0.1*range_x,max(final_time)+0.1*range_x])
    #axes.set_ylim([min(final_h)-0.1*range_y,max(final_h)+0.1*range_y])
    plt.xlabel('time')
    plt.ylabel('h')
    plt.title('height')
    plt.grid(True)
    plt.savefig("height.png")
    #plt.show()

def plot_xy(final_x,final_y):
    plt.plot(final_x, final_y)
    range_x = max(final_x)-min(final_x)
    range_y = max(final_y)-min(final_y)
    axes = plt.gca()
    axes.set_xlim([min(final_x)-0.1*range_x,max(final_x)+0.1*range_x])
    axes.set_ylim([min(final_y)-0.1*range_y,max(final_y)+0.1*range_y])
    #plt.gca().invert_yaxis()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('x y')
    plt.grid(True)
    plt.savefig("xy.png")
    #plt.show()
