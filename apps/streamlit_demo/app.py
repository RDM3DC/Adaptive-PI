import streamlit as st
from pi_a.core import triangle_angle_sum_pia, is_pia_cyclic
from pi_a.models import GaussianBump

st.title("Adaptive π (πₐ) — Live Demo")

Ax = st.slider("Ax", -1.0, 2.0, 0.0, 0.01)
Ay = st.slider("Ay", -1.0, 2.0, 0.0, 0.01)
Bx = st.slider("Bx", -1.0, 2.0, 1.0, 0.01)
By = st.slider("By", -1.0, 2.0, 0.0, 0.01)
Cx = st.slider("Cx", -1.0, 2.0, 0.2, 0.01)
Cy = st.slider("Cy", -1.0, 2.0, 0.8, 0.01)

K0 = st.slider("Gaussian K0", -0.1, 0.1, 0.02, 0.001)
x0 = st.slider("x0", -1.0, 2.0, 0.5, 0.01)
y0 = st.slider("y0", -1.0, 2.0, 0.35, 0.01)
sig = st.slider("sigma", 0.05, 1.0, 0.3, 0.01)

A,B,C = (Ax,Ay), (Bx,By), (Cx,Cy)
Kb = GaussianBump(K0=K0, x0=x0, y0=y0, sigma=sig)
S = triangle_angle_sum_pia(A,B,C, Kb)

st.write("**πₐ triangle angle sum:**", S)

Xx = st.slider("Xx", -1.0, 2.0, 0.9, 0.01)
Xy = st.slider("Xy", -1.0, 2.0, 0.2, 0.01)
X = (Xx, Xy)

st.write("**πₐ-cyclic(A,B,C,X)?**", is_pia_cyclic(A,B,C, X, Kb))
