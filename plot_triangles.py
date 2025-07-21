# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Arc#, Circle
import os
os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"  # Force software rendering
import CoolProp.CoolProp as CP


from CoolProp import AbstractState


# Define the backend and fluid
fluid = 'Air'
AS = AbstractState("HEOS", fluid)  
# AS = AbstractState("TTSE", "Air")  

def CP_AS_PT_HS(AS,P,T):
    AS.update(CP.PT_INPUTS, P, T)
    return AS.hmass(),AS.smass()

def CP_AS_HS_P(AS,H,S):
    AS.update(CP.HmassSmass_INPUTS, H,S,)
    return AS.p()

def CP_AS_HP_S(AS,H,P):
    AS.update(CP.HmassP_INPUTS, H,P)
    return AS.smass()






# todo
# nice angles in plot
# do this in jupyter so that everybody can open it on the phone?
# gitlab, QR code

import plotting_functions as pf

def velocity_triangle(u,FlowCoef,WorkCoef,Reaction,AS):
    """calculates all relevant values of the velocity triangles"""
    u1 = u  # Assumption: Constant meanline radius (axial machine)
    u2 = u
    cm1 = FlowCoef * u1 # FlowCoef = cm1 / u 
    # cx1 = cm1 # Assumption: no radial component at inlet (axial machine)
    cm = cm1 # constant axial velocity
    
    cu1 = u * (1 - Reaction - 0.5*WorkCoef)
    cu2 = u * (1 - Reaction + 0.5*WorkCoef)

    wu1 = cu1 - u1 
    wu2 = cu2 - u2
    
    alpha1 = np.degrees(np.arctan2(cu1,cm))
    beta1 = np.degrees(np.arctan2(wu1,cm))

    alpha2 = np.degrees(np.arctan2(cu2,cm))
    beta2 = np.degrees(np.arctan2(wu2,cm))

    c1 = (cm**2+cu1**2)**0.5
    w1 = (cm**2+wu1**2)**0.5
    c2 = (cm**2+cu2**2)**0.5
    w2 = (cm**2+wu2**2)**0.5



    delta_h12t = WorkCoef * u**2
    delta_h12 = delta_h12t  - c2**2/2 + c1**2/2
    if Reaction != 0:
        delta_h_stage = delta_h12 / Reaction
    else:
        delta_h_stage = delta_h12
        
    pt_in = 30e5
    Tt_in = 1800

    eta_poly = 0.8
    
    def sout_from_eta_poly(WorkCoef, eta,s_in,s_in_h):
        if WorkCoef < 0: #Turbine 
            # eta_poly = (s_in-s_in_h) / (s_out-s_in_h)     # turbine
            s_out = (s_in-s_in_h) / eta  + s_in_h
        elif WorkCoef >= 0: #Compressor 
            # eta_poly = (s_in_h-s_out) / (s_in_h-s_in)     # compressor
            s_out = s_in_h - (eta_poly * (s_in_h-s_in)) 
        return s_out
    
    if WorkCoef <= 0: #Turbine
        # first stator(0,1) -> then rotor (1,2) 
        p0t = pt_in
        T0t = Tt_in
        
        # h0t = CP.PropsSI('HMASS', 'P', p0t, 'T', T0t, fluid)  
        # s0 = CP.PropsSI('SMASS', 'P', p0t, 'T', T0t, fluid) 
        h0t,s0 = CP_AS_PT_HS(AS,p0t,T0t)
        c0 = c2
        h0 = h0t - c0**2/2
        # p0 = CP.PropsSI('P','HMASS',  h0, 'SMASS', s0, fluid)
        p0 = CP_AS_HS_P(AS,h0,s0)
        
        
        # # delta_h_stage = h0-h2
        # # delta_h12 = h1-h2 -> h2 = h1 - delta_h12
        # # delta_h_stage = h0- (h1 - delta_h12) = h0 - h1 + delta_h12
        # h1 = h0 + delta_h12 - delta_h_stage
        
        # r = h2-h1 / h2 - h0
        # delta_h_stage = h2-h0
        # delta_h12 = h2-h1 -> h2 = h1 + delta_h12
        # delta_h_stage = (h1 + delta_h12) -h0 =
        h1 = delta_h_stage + h0 - delta_h12
        h1t = h1 + c1**2/2
        
        # s0h = CP.PropsSI('SMASS', 'P', p0, 'HMASS', h1, fluid)
        s0h = CP_AS_HP_S(AS,h1,p0)
        s1 = sout_from_eta_poly(WorkCoef,eta_poly,s0,s0h)
        
    elif WorkCoef > 0: #Compressor
        # first rotor(1,2) -> then stator (2,3) 
    
        p1t = pt_in
        T1t = Tt_in
        
        # h1t = CP.PropsSI('HMASS', 'P', p1t, 'T', T1t, fluid)  
        # s1 = CP.PropsSI('SMASS', 'P', p1t, 'T', T1t, fluid) 
        h1t,s1 = CP_AS_PT_HS(AS,p1t,T1t)
        
    
    h1 = h1t - c1**2/2

    # p1 = CP.PropsSI('P','HMASS',  h1, 'SMASS', s1, fluid)
    p1 = CP_AS_HS_P(AS,h1,s1)
    
    
    # p1t = CP.PropsSI('P','HMASS',  h1t, 'SMASS', s1, fluid)
    h2t = delta_h12t + h1t
    h2 = delta_h12 + h1
    
    # s1h = CP.PropsSI('SMASS', 'P', p1, 'HMASS', h2, fluid)
    s1h = CP_AS_HP_S(AS,h2,p1)
    s2 = sout_from_eta_poly(WorkCoef,eta_poly,s1,s1h)
    # p2 = CP.PropsSI('P','HMASS',  h2, 'SMASS', s2, fluid)
    p2 = CP_AS_HS_P(AS,h2,s2)
    h2t = h2 + c2**2/2
    # p2t = CP.PropsSI('P','HMASS',  h2t, 'SMASS', s2, fluid)
    
    if WorkCoef > 0: #Compressor

        # Reaction = delta_h12/(h3-h1)
        h3 = delta_h12/Reaction + h1
        c3 = c1
        h3t = h3 + c3**2/2
        
        # s2h = CP.PropsSI('SMASS', 'P', p2, 'HMASS', h3, fluid)
        s2h = CP_AS_HP_S(AS,h3,p2)
        s3 = sout_from_eta_poly(WorkCoef,eta_poly,s2,s2h)
        # p3 = CP.PropsSI('P','HMASS',  h3, 'SMASS', s3, fluid)
        p3 = CP_AS_HS_P(AS,h3,s3)
        # for easier variable handling: 
        h0=h3
        h0t=h3t
        s0=s3
        p0=p3
        

    hs_dict = {
                # 'c1':c1,'c2':c2,'w1':w1,'w2':w2,'cm':cm,'cu1':cu1,'cu2':cu2,'wu1':wu1,'wu2':wu2,
               # 'alpha1':alpha1,'beta1':beta1,'alpha2':alpha2,'beta2':beta2,
               # 'p0t':p0t,'p0':p0,'p1t':p1t,'p1':p1,'p2t':p2t,'p2':p2,
               'p0':p0,'p1':p1,'p2':p2,
               'h0t':h0t,'h0':h0,'h1t':h1t,'h1':h1,'h2t':h2t,'h2':h2,
               's0':s0,'s1':s1,'s2':s2}



    # return c1,c2,w1,w2,cm,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2,p0t,p0,p1t,p1,p2t,p2
    return c1,c2,w1,w2,cm,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2, hs_dict
    # return res_dict











fig = plt.figure()
# fig.suptitle("Velocity Triangle")


gs = GridSpec(8, 11)


ax0 = fig.add_subplot(gs[0, :2])
ax1 = fig.add_subplot(gs[1, :2])
ax2 = fig.add_subplot(gs[2, :2])
ax3 = fig.add_subplot(gs[3, :2])
ax4 = fig.add_subplot(gs[4:, :2])
ax5 = fig.add_subplot(gs[:4, 3:7]) # triangle
ax6 = fig.add_subplot(gs[:4, 7:11]) # hs
ax7 = fig.add_subplot(gs[4:, 3:]) # camberlines



# Add sliders to modify plots
slider_u = Slider(ax0, 'u in m/s', 100, 300, valinit=200, valstep=10)
slider_flow = Slider(ax1, 'FlowCoef $\\frac{{c_x}}{{u}}$', 0.3, 1, valinit=0.5, valstep=0.05)
slider_work = Slider(ax2, 'WorkCoef $\\frac{{c_{{u2}} - c_{{u1}}}}{{u}}$', -2.5, 0.9, valinit=-1.2, valstep=0.05)
slider_reaction = Slider(ax3, 'Reaction', -0.1, 1.0, valinit=0.5, valstep=0.01)

u=200
c1,c2,w1,w2,cm,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2, hs_dict = velocity_triangle(u,0.5,0.5,0.5,AS)    



### Velocity triangles
ax5.set_title('velocity triangles')
ax5.axis('equal')
ax5.set_axis_off()

l_u = pf.draw_arrow(ax5,[u,0],u,0,'u','C0',linestyle='-')
l_c1 = pf.draw_arrow(ax5,[u,0],cu1,cm,'c1','C1',linestyle='-')
l_w1 = pf.draw_arrow(ax5,[0,0],wu1,cm,'w1','C2',linestyle='-',horizontalalignment='right')
l_c2 = pf.draw_arrow(ax5,[u,0],cu2,cm,'c2','C1',linestyle='--')
l_w2 = pf.draw_arrow(ax5,[0,0],wu2,cm,'w2','C2',linestyle='--',horizontalalignment='right')






### hs diagram
ax6.set_title('enthalpy - entropy diagram')
ax6.set_xticklabels([])  # hides x-axis labels
ax6.set_yticklabels([])  # hides y-axis labels
ax6.set_xlabel('Entropy')
ax6.set_ylabel('Enthalpy')
ax6.plot()


hs_line = ax6.plot([hs_dict['s0'],hs_dict['s1'],hs_dict['s2']],[hs_dict['h0'],hs_dict['h1'],hs_dict['h2']],marker='s',linestyle='-') 
h0t_line = ax6.plot([hs_dict['s0'],hs_dict['s0']],[hs_dict['h0'],hs_dict['h0t']],marker='s',linestyle='-') 
h1t_line = ax6.plot([hs_dict['s1'],hs_dict['s1']],[hs_dict['h1'],hs_dict['h1t']],marker='s',linestyle='-') 
h2t_line = ax6.plot([hs_dict['s2'],hs_dict['s2']],[hs_dict['h2'],hs_dict['h2t']],marker='s',linestyle='-') 

hmin = min(hs_dict['h2'],hs_dict['h1'],)
hmax = max(hs_dict['h2t'],hs_dict['h1t'])
h_delta = hmax - hmin
range_fraction = 0.2
h_range = np.linspace(hmin-h_delta*range_fraction,hmax+h_delta*range_fraction,101)

vectorized_func = np.vectorize(lambda H: CP_AS_HP_S(AS, H, hs_dict['p0']))
p0_isoline_data = vectorized_func(h_range)
vectorized_func = np.vectorize(lambda H: CP_AS_HP_S(AS, H, hs_dict['p1']))
p2_isoline_data = vectorized_func(h_range)
vectorized_func = np.vectorize(lambda H: CP_AS_HP_S(AS, H, hs_dict['p2']))
p1_isoline_data = vectorized_func(h_range)
p0_isoline = ax6.plot(p0_isoline_data,h_range,color='k',alpha=0.5)
p1_isoline = ax6.plot(p1_isoline_data,h_range, color='k',alpha=0.5)
p2_isoline = ax6.plot(p2_isoline_data,h_range, color='k',alpha=0.5)

smin = min(hs_dict['s0'],hs_dict['s1'])
smax = max(hs_dict['s0'],hs_dict['s2'])
s_delta = smax - smin
s_axlim_min = smin -s_delta*range_fraction
s_axlim_max = smax +s_delta*range_fraction
    
ax6.set_ylim(h_range.min(),h_range.max())
ax6.set_xlim(s_axlim_min,s_axlim_max)

h1t_text = ax6.text(hs_dict['s1']+s_delta*0.01,hs_dict['h1t']+h_delta*0.01,'h1t', size='large')
h1_text = ax6.text(hs_dict['s1']+s_delta*0.01,hs_dict['h1']+h_delta*0.01,'h1', size='large')
h2t_text = ax6.text(hs_dict['s2']+s_delta*0.01,hs_dict['h2t']+h_delta*0.01,'h2t', size='large')
h2_text = ax6.text(hs_dict['s2']+s_delta*0.01,hs_dict['h2']+h_delta*0.01,'h2', size='large')



### Camber line arcs
# fixed parameters
ax7.set_title('camber lines')
ax7.axis('equal')
ax7.set_axis_off()

# initial plot 
l_u_camber = pf.draw_arrow(ax7,[u,0],u,0,'u','C0',linestyle='-')
l_c1_camber = pf.draw_arrow(ax7,[u,0],cu1,cm,'c1','C1',linestyle='-')
l_w1_camber = pf.draw_arrow(ax7,[0,0],wu1,cm,'w1','C2',linestyle='-',horizontalalignment='right')

P0     = (0, 0)
h      = -100      # vertical rise
# compute initial geometry
P1, C, R, t1, t2, is_line = pf.compute_arc_geometry(P0, beta1, beta2, h)

# create both patches up front
arc1_patch  = Arc(xy=P0, width=0, height=0, angle=0, theta1=0, theta2=0, edgecolor='k', lw=5)
ax7.add_patch(arc1_patch)

# Create a second arc patch for the offset arc
arc2_patch = Arc(xy=P0, width=0, height=0, angle=0,theta1=0, theta2=0, edgecolor='k', lw=5, ls='-')
line_patch, = ax7.plot([], [], '-', color='k', lw=5)
ax7.add_patch(arc2_patch)

circle_patches = []
circle_patches = pf.draw_camber_circles(fig,ax7,is_line,arc1_patch,arc2_patch,line_patch,P0,P1,C,R,t1,t2,h,circle_patches)

l_u2_camber = pf.draw_arrow(ax7,[P1[0]+cu2,P1[1]-cm],u,0,'u','C0',linestyle='-')
l_c2_camber = pf.draw_arrow(ax7,[P1[0]+cu2,P1[1]-cm],cu2,cm,'c2','C1',linestyle='--')
l_w2_camber = pf.draw_arrow(ax7,[P1[0]+wu2,P1[1]-cm],wu2,cm,'w2','C2',linestyle='--',horizontalalignment='right')





### text field
text_ax4 = ax4.text(0.5, 0.5, pf.angle_text(u,c1,c2,w1,w2,cm,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2), size='large', va="center", ha="center")
ax4.set_axis_off()


def update_plot(val):
    """ function which updates the plot"""
    u = slider_u.val

    FlowCoef = slider_flow.val
    WorkCoef = slider_work.val
    Reaction = slider_reaction.val
    
    c1,c2,w1,w2,cm,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2, hs_dict = velocity_triangle(u,FlowCoef,WorkCoef,Reaction,AS)    
    
    ### velocity triangles
    pf.redraw_arrow(l_u,[u,0],u,0)
    pf.redraw_arrow(l_c1,[u,0],cu1,cm)
    pf.redraw_arrow(l_w1,[0,0],wu1,cm)
    pf.redraw_arrow(l_c2,[u,0],cu2,cm)
    pf.redraw_arrow(l_w2,[0,0],wu2,cm)
    
    
    ### text field
    text_ax4.set_text(pf.angle_text(u,c1,c2,w1,w2,cm,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2))



    ### hs diagram
    
    
    # hs_line = ax6.plot([hs_dict['s0'],hs_dict['s1'],hs_dict['s2']],[hs_dict['h0'],hs_dict['h1'],hs_dict['h2']],marker='s',linestyle='-') 
    
    hs_line[0].set_xdata([hs_dict['s0'],hs_dict['s1'],hs_dict['s2']])
    hs_line[0].set_ydata([hs_dict['h0'],hs_dict['h1'],hs_dict['h2']])
    
    h0t_line[0].set_xdata([hs_dict['s0'],hs_dict['s0']])
    h0t_line[0].set_ydata([hs_dict['h0'],hs_dict['h0t']])
    h1t_line[0].set_xdata([hs_dict['s1'],hs_dict['s1']])
    h1t_line[0].set_ydata([hs_dict['h1'],hs_dict['h1t']])
    h2t_line[0].set_xdata([hs_dict['s2'],hs_dict['s2']])
    h2t_line[0].set_ydata([hs_dict['h2'],hs_dict['h2t']])
    
    
    
    
    smin = min(hs_dict['s0'],hs_dict['s1'])
    smax = max(hs_dict['s0'],hs_dict['s2'])
    s_delta = smax - smin
    s_axlim_min = smin -s_delta*range_fraction
    s_axlim_max = smax +s_delta*range_fraction
    # s_range = np.linspace(smin-s_delta*s_fraction,smax+s_delta*s_fraction,101)
    # p0_isoline_data = CP.PropsSI('HMASS', 'P', hs_dict['p0'], 'SMASS', s_range, 'Air')
    # p1_isoline_data = CP.PropsSI('HMASS', 'P', hs_dict['p1'], 'SMASS', s_range, 'Air')
    # p2_isoline_data = CP.PropsSI('HMASS', 'P', hs_dict['p2'], 'SMASS', s_range, 'Air')
    
    hmin = min(hs_dict['h2'],hs_dict['h1'],)
    hmax = max(hs_dict['h2t'],hs_dict['h1t'])
    h_delta = hmax - hmin
    
    h_range = np.linspace(hmin-h_delta*range_fraction,hmax+h_delta*range_fraction,101)

    vectorized_func = np.vectorize(lambda H: CP_AS_HP_S(AS, H, hs_dict['p0']))
    p0_isoline_data = vectorized_func(h_range)
    vectorized_func = np.vectorize(lambda H: CP_AS_HP_S(AS, H, hs_dict['p1']))
    p1_isoline_data = vectorized_func(h_range)
    vectorized_func = np.vectorize(lambda H: CP_AS_HP_S(AS, H, hs_dict['p2']))
    p2_isoline_data = vectorized_func(h_range)
   
    
    
    p0_isoline[0].set_ydata(h_range)
    p0_isoline[0].set_xdata(p0_isoline_data)
    p1_isoline[0].set_ydata(h_range)
    p1_isoline[0].set_xdata(p1_isoline_data)
    p2_isoline[0].set_ydata(h_range)
    p2_isoline[0].set_xdata(p2_isoline_data)
    ax6.set_ylim(h_range.min(),h_range.max())
    ax6.set_xlim(s_axlim_min,s_axlim_max)


    h1t_text.set_x(hs_dict['s1']+s_delta*0.01)
    h1t_text.set_y(hs_dict['h1t']+h_delta*0.01)
    h1_text.set_x(hs_dict['s1']+s_delta*0.01)
    h1_text.set_y(hs_dict['h1']+h_delta*0.01)
    h2t_text.set_x(hs_dict['s2']+s_delta*0.01)
    h2t_text.set_y(hs_dict['h2t']+h_delta*0.01)
    h2_text.set_x(hs_dict['s2']+s_delta*0.01)
    h2_text.set_y(hs_dict['h2']+h_delta*0.01)
    
  
    
    
    
    
    ### camber
    pf.redraw_arrow(l_u_camber,[u,0],u,0)
    pf.redraw_arrow(l_c1_camber,[u,0],cu1,cm)
    pf.redraw_arrow(l_w1_camber,[0,0],wu1,cm)
    
    global circle_patches
    P1, C, R, t1, t2, is_line = pf.compute_arc_geometry(P0, beta1, beta2, h)
    pf.draw_camber_circles(fig,ax7,is_line,arc1_patch,arc2_patch,line_patch,P0,P1,C,R,t1,t2,h,circle_patches)
    
    pf.redraw_arrow(l_u2_camber,[P1[0]+cu2,P1[1]-cm],u,0)
    pf.redraw_arrow(l_c2_camber,[P1[0]+cu2,P1[1]-cm],cu2,cm)
    pf.redraw_arrow(l_w2_camber,[P1[0]+wu2,P1[1]-cm],wu2,cm)
    return 


slider_u.on_changed(update_plot)
slider_flow.on_changed(update_plot)
slider_work.on_changed(update_plot)
slider_reaction.on_changed(update_plot)


plt.show()