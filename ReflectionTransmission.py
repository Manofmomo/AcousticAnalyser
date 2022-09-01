#!/usr/bin/env python
# coding: utf-8

# In[16]:


import sympy
from sympy import I
rho1, Area1, E1, I1 = sympy.symbols("rho1, Area1, E1, I1")
rho2, Area2, E2, I2 = sympy.symbols("rho2, Area2, E2, I2")

x = sympy.symbols("x")
a, b, c, d= sympy.symbols("a, b, c, d")
ad, bd, cd, dd= sympy.symbols("ad, bd, cd, dd")
ks= sympy.symbols("Ks")
w = sympy.symbols("w")

k1 = (rho1*Area1/(E1*I1))**(0.25)*(w**0.5)
k2 = (rho2*Area2/(E2*I2))**(0.25)*(w**0.5)


# In[17]:


def subs_nums(eq):
    return eq.subs(ks,1).subs(rho1,1).subs(Area1,1).subs(E1,1).subs(I1,1).subs(rho2,1).subs(Area2,1).subs(E2,2).subs(I2,2).expand()


# In[18]:


incoming_wave = sympy.exp(-I*k1*x)


# In[19]:


v1 = incoming_wave + c*sympy.exp(I*k1*x) + d*sympy.exp(k1*x)
v2 = a*sympy.exp(-I*k2*x) + b*sympy.exp(-k2*x)


# In[20]:


v = sympy.Piecewise((v1,x<0),(v2,x>0))


# In[21]:


vdiffs= [None]*4
vdiffs[0] = v
for i in range(3):
    vdiffs[i+1] = sympy.diff(vdiffs[i],x)


# In[22]:


eq1 = sympy.limit(vdiffs[0], x, 0, dir='-')-sympy.limit(vdiffs[0], x, 0, dir='+')
eq2 = sympy.limit(vdiffs[1], x, 0, dir='-')-sympy.limit(vdiffs[1], x, 0, dir='+')
eq3 = E1*I1*sympy.limit(vdiffs[2], x, 0, dir='-')-E2*I2*sympy.limit(vdiffs[2], x, 0, dir='+')
eq4 = E1*I1*sympy.limit(vdiffs[3], x, 0, dir='-')+ks*sympy.limit(vdiffs[0], x, 0, dir='+')-E2*I2*sympy.limit(vdiffs[3], x, 0, dir='+')


# In[23]:





# In[24]:


soln = sympy.solve((eq1,eq2,eq3,eq4),(a,b,c,d))


# In[25]:


a_sol = soln[a].expand()
b_sol = soln[b].expand()
c_sol = soln[c].expand()
d_sol = soln[d].expand()


# In[ ]:





# In[ ]:





# In[ ]:





# In[26]:


incoming_wave = sympy.exp(-1*k1*x)


# In[27]:


v1 = incoming_wave + cd*sympy.exp(I*k1*x) + dd*sympy.exp(k1*x)
v2 = ad*sympy.exp(-I*k2*x) + bd*sympy.exp(-k2*x)


# In[28]:


v = sympy.Piecewise((v1,x<0),(v2,x>0))


# In[29]:


vdiffs= [None]*4
vdiffs[0] = v
for i in range(3):
    vdiffs[i+1] = sympy.diff(vdiffs[i],x)


# In[30]:


eq1 = sympy.limit(vdiffs[0], x, 0, dir='-')-sympy.limit(vdiffs[0], x, 0, dir='+')
eq2 = sympy.limit(vdiffs[1], x, 0, dir='-')-sympy.limit(vdiffs[1], x, 0, dir='+')
eq3 = E1*I1*sympy.limit(vdiffs[2], x, 0, dir='-')-E2*I2*sympy.limit(vdiffs[2], x, 0, dir='+')
eq4 = E1*I1*sympy.limit(vdiffs[3], x, 0, dir='-')+ks*sympy.limit(vdiffs[0], x, 0, dir='+')-E2*I2*sympy.limit(vdiffs[3], x, 0, dir='+')


# In[ ]:





# In[32]:


soln_2 = sympy.solve((eq1,eq2,eq3,eq4),(ad,bd,cd,dd))


# In[ ]:





# In[33]:


ad_sol = soln_2[ad].expand()
bd_sol = soln_2[bd].expand()
cd_sol = soln_2[cd].expand()
dd_sol = soln_2[dd].expand()


# In[ ]:





# In[34]:


transmission_matrix = sympy.Matrix([[a_sol,b_sol],[ad_sol, bd_sol]])
reflection_matrix = sympy.Matrix([[c_sol,d_sol],[cd_sol, dd_sol]])


# In[ ]:





# In[ ]:





# In[37]:


import pickle 
filehandler = open('reflection_transmission_matrix', 'wb') 
pickle.dump((reflection_matrix,transmission_matrix), filehandler)
filehandler.close()


# In[ ]:





# In[46]:


print(subs_nums(reflection_matrix).simplify())


# In[44]:


print(subs_nums(transmission_matrix).simplify())

