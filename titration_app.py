import streamlit as st
import numpy as np
from scipy.optimize import fsolve
import plotly.express as px

# Calculate affect of weak acid on pH
# Done by simultaneously solving K_a and K_w
# x is concentration of HA hydrolyzed
# y is concentration of H_2O self-hydrolyzed
@st.cache_data(max_entries=2000)
def solve_weak_acid(K_a, c_h, c_oh, c_HA, c_A):
    guess=[np.sqrt(K_a*c_HA), -np.sqrt(K_a*c_HA)/100]
    def equations(vars):
        x = vars[0]
        y = vars[1]

        eq1 = x**2 + x*(c_A+K_a+c_h) + (c_h*c_A-c_HA*K_a) + y*c_A + x*y
        eq2 = c_oh*x + (c_oh+c_h)*y + x*y + y**2
        return [eq1, eq2]

    # Solve the system numerically
    return fsolve(equations, guess)

def calc_pH(K_a, c0_acid, V0_acid, c0_base, V0_base):
    # Convert ml to l
    V0_acid_l = np.array(V0_acid/1000)
    V0_base_l = V0_base/1000

    # Final volume
    V1 = V0_acid_l+V0_base_l

    # Calculate initial number of moles
    n0_acid = c0_acid*V0_acid_l*np.ones_like(V1)
    n0_base = c0_base*V0_base_l*np.ones_like(V1)
    
    # Calculate number of moles after acid-base reaction
    # Depends on if all base is used up or not
    no_excess_base = n0_acid>=n0_base
    n_HA = np.where(no_excess_base, n0_acid-n0_base, 0)
    n_A = np.where(no_excess_base, n0_base, n0_acid)
    n_base=np.where(no_excess_base, 0, n0_base-n0_acid)
         
    # Calculate concentrations after acid-base reaction
    c_HA = n_HA/V1
    c_A = n_A/V1
    c_base = n_base/V1

    # If there's excess base, calculate affect on pH
    c_oh = 1e-7+c_base
    c_h = 1e-14/c_oh

    # Calculate affect of weak acid on pH
    # Done by simultaneously solving k_A and k_w
    # x is concentration of H_3O+ gained from hydrolysis of HA
    # y is concentration of H_3O+ gained from self-hydrolysis of water
    x = np.zeros_like(V1)
    y = np.zeros_like(V1)

    if V1.size==1:
        x, y = solve_weak_acid(K_a, c_h, c_oh, c_HA, c_A)
    else:
        for i in np.arange(V1.size):
            x[i], y[i] = solve_weak_acid(K_a, c_h[i], c_oh[i], c_HA[i], c_A[i])

    # Add change in H_3O+ and return pH
    return -np.log10(c_h+x+y)

def plot_pH(K_a, c0_acid, V0_acid, c0_base, annotate=True):
    # Double equivalence
    V_max_base = 2*c0_acid*V0_acid/c0_base
    V0_base = np.arange(V_max_base+1)


    pH = calc_pH(K_a, c0_acid, V0_acid, c0_base, V0_base)

    fig = px.line(x=V0_base, y=pH, 
        labels={"x":"Volume of strong base (ml)", "y": "pH"},
        range_y=[2,14],
        template="seaborn")
    
    fig.add_hline(y=7, line_color="magenta", line_dash="dash", annotation_text="    Neutral pH", annotation_position="top left", opacity=1, line_width=2, annotation_font_size=20)

    fig.update_layout(
        xaxis_title_font_size=20,
        yaxis_title_font_size=20,
        xaxis_tickfont_size=16,
        yaxis_tickfont_size=16,
        hoverlabel_font_size=16,
        xaxis=dict(fixedrange=True),
        yaxis=dict(fixedrange=True)
        )

    return fig

with st.sidebar:
    with st.form("params"):
        K_a = st.number_input("Weak acid dissociation constant K_a", value=1.8e-5, min_value=1e-14, format="%.2e", step=1e-6)

        c0_acid = st.number_input("Weak acid concentration (M)", value=0.1, min_value=1e-7)
        V0_acid = st.number_input("Weak acid volume (ml)", value=40, min_value=1)

        c0_base = st.number_input("Strong base concentration (M)", value=0.1, min_value=1e-7)
        V_max = 2*c0_acid*V0_acid/c0_base

        run = st.form_submit_button("Run")

fig = plot_pH(K_a, c0_acid, V0_acid, c0_base, annotate=False)
plot = st.plotly_chart(fig, theme=None)

# propanoic
# K_a = 1.3e-5

# acetic
# K_a = 1.8e-5

# nitrous
# K_a = 4e-4

# hydrafluouric
# K_a = 6.3e-4
