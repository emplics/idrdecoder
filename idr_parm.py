from tensorflow.python.keras import backend as K
#import tensorflow as tf
import numpy as np

para1 = 800
para2 = 21
encoding_dim = 20
norm1 = 1   
norm3 = 2   
norm4 = 0.05
weight     = 100 
weightpg   = 400000 
threshold  = 0.70 
threshold_pg = 0.55 
diff_range = 3  
diff_rate  = 0.005 
diff_rate3 = 0.000
numpg = 87 
pc1 = 1
pc2 = 2
scorescale = 1
showno = 1 #1:seq no, 2:seq name
#
def pearson_r(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x, axis=0)
    my = K.mean(y, axis=0)
    xm, ym = x - mx, y - my
    r_num = K.sum(xm * ym)
    x_square_sum = K.sum(xm * xm)
    y_square_sum = K.sum(ym * ym)
    r_den = K.sqrt(x_square_sum * y_square_sum)
    r = 0.
    if r_den != 0.:
      r = r_num / r_den
    return K.mean(r)

def pearson_r_loss(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x, axis=0)
    my = K.mean(y, axis=0)
    xm, ym = x - mx, y - my
    r_num = K.sum(xm * ym)
    x_square_sum = K.sum(xm * xm)
    y_square_sum = K.sum(ym * ym)
    r_den = K.sqrt(x_square_sum * y_square_sum)
    r = 0.
    if r_den != 0.:
      r = r_num / r_den
    return 1 - K.square(K.mean(r))
