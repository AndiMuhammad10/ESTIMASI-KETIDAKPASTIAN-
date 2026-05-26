import streamlit as st
import math
import pandas as pd
import re

# ═══════════════════════════════════════════
# STREAMLIT CONFIG
# ═══════════════════════════════════════════
st.set_page_config(
    page_title="Estimasi Ketidakpastian Titrimetri",
    page_icon="⚗️",
    layout="wide"
)

# ═══════════════════════════════════════════
# CUSTOM CSS
# ═══════════════════════════════════════════
st.markdown("""
<style>

/* Background utama */
.stApp {
    background: linear-gradient(
        135deg,
        #14001f 0%,
        #2b0a3d 25%,
        #4c1d95 60%,
        #7c3aed 100%
    );
    color: #f5e9ff;
}

/* Semua teks */
html, body, [class*="css"] {
    color: #f5e9ff !important;
    font-family: 'Segoe UI', sans-serif;
}

/* Judul */
h1, h2, h3, h4, h5, h6 {
    color: white !important;
    text-shadow: 0 0 10px rgba(255,255,255,0.25);
}

/* Sidebar */
section[data-testid="stSidebar"] {
    background: linear-gradient(
        180deg,
        #1e0033 0%,
        #3c096c 100%
    );
}

/* Input */
.stTextInput input,
.stNumberInput input,
textarea {
    background-color: rgba(255,255,255,0.08) !important;
    color: white !important;
    border: 1px solid #a855f7 !important;
    border-radius: 12px !important;
}

/* Selectbox */
div[data-baseweb="select"] > div {
    background-color: rgba(255,255,255,0.08) !important;
    color: white !important;
    border: 1px solid #a855f7 !important;
    border-radius: 12px !important;
}

/* Tombol */
.stButton button {
    background: linear-gradient(
        90deg,
        #7c3aed,
        #a855f7
    ) !important;

    color: white !important;
    border: none !important;
    border-radius: 14px !important;
    font-weight: bold;
}

/* Metric */
[data-testid="metric-container"] {
    background: rgba(255,255,255,0.08);
    border-radius: 16px;
    padding: 15px;
}

/* Latex */
.katex {
    color: white !important;
}

/* Dataframe */
[data-testid="stDataFrame"] {
    background: rgba(255,255,255,0.05);
    border-radius: 14px;
}

</style>
""", unsafe_allow_html=True)

# ═══════════════════════════════════════════
# JUDUL APP
# ═══════════════════════════════════════════
st.title("⚗️ Estimasi Ketidakpastian Titrimetri")
st.caption("Kalkulasi µ komponen dan ketidakpastian gabungan U₉₅ secara otomatis")

# ─────────────────────────────────────────────
# DATABASE TABEL PERIODIK
# ─────────────────────────────────────────────
PERIODIC_DB = {
    "H":  {"ar": 1.00794,   "qu": 7e-5},
    "He": {"ar": 4.002602,  "qu": 2e-6},
    "Li": {"ar": 6.941,     "qu": 2e-3},
    "B":  {"ar": 10.811,    "qu": 7e-5},
    "C":  {"ar": 12.0107,   "qu": 8e-5},
    "N":  {"ar": 14.0067,   "qu": 1e-4},
    "O":  {"ar": 15.9994,   "qu": 3e-5},
    "F":  {"ar": 18.9984,   "qu": 1e-6},
    "Na": {"ar": 22.98977,  "qu": 2e-5},
    "Mg": {"ar": 24.3050,   "qu": 6e-4},
    "Al": {"ar": 26.98154,  "qu": 8e-6},
    "Si": {"ar": 28.0855,   "qu": 3e-4},
    "P":  {"ar": 30.97376,  "qu": 2e-6},
    "S":  {"ar": 32.065,    "qu": 5e-4},
    "Cl": {"ar": 35.453,    "qu": 2e-3},
    "K":  {"ar": 39.0983,   "qu": 1e-4},
    "Ca": {"ar": 40.078,    "qu": 4e-4},
    "Mn": {"ar": 54.9380,   "qu": 1e-5},
    "Fe": {"ar": 55.845,    "qu": 2e-3},
    "Co": {"ar": 58.9332,   "qu": 2e-5},
    "Ni": {"ar": 58.6934,   "qu": 2e-4},
    "Cu": {"ar": 63.546,    "qu": 3e-4},
    "Zn": {"ar": 65.38,     "qu": 2e-3},
    "Br": {"ar": 79.904,    "qu": 1e-3},
    "Ag": {"ar": 107.8682,  "qu": 2e-4},
    "I":  {"ar": 126.9045,  "qu": 3e-4},
}

# ─────────────────────────────────────────────
# PARSER FORMULA
# ─────────────────────────────────────────────
def parse_formula(formula: str):

    result = {}

    def parse_group(s, mult=1):

        counts = {}

        pattern = re.compile(
            r'([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)'
        )

        stack = [{}]

        for m in pattern.finditer(s):

            el, n_str, open_p, close_p, close_n = m.groups()

            if el:

                n = int(n_str) if n_str else 1

                stack[-1][el] = stack[-1].get(el, 0) + n

            elif open_p:

                stack.append({})

            elif close_p:

                n = int(close_n) if close_n else 1

                top = stack.pop()

                for k, v in top.items():
                    stack[-1][k] = stack[-1].get(k, 0) + v * n

        for k, v in stack[0].items():
            counts[k] = counts.get(k, 0) + v * mult

        return counts

    parts = formula.replace("·", ".").split(".")

    for part in parts:

        m = re.match(r'^(\d+)?(.*)', part.strip())

        mult = int(m.group(1)) if m.group(1) else 1

        sub = m.group(2).strip()

        if sub:

            sub_counts = parse_group(sub, mult)

            for k, v in sub_counts.items():
                result[k] = result.get(k, 0) + v

    return result

# ─────────────────────────────────────────────
# HITUNG µ_BM
# ─────────────────────────────────────────────
def calc_mu_BM(formula, n_equiv=1):

    parsed = parse_formula(formula)

    sqrt3 = math.sqrt(3)

    bm = 0
    sum_uc2 = 0

    rows = []

    unknown = [el for el in parsed if el not in PERIODIC_DB]

    if unknown:
        return None, None, None, None, None, f"Unsur tidak dikenali: {', '.join(unknown)}"

    for el, n in parsed.items():

        ar = PERIODIC_DB[el]["ar"]
        qu = PERIODIC_DB[el]["qu"]

        n_ar = n * ar

        muc = n * qu / sqrt3

        muc2 = muc ** 2

        bm += n_ar
        sum_uc2 += muc2

        rows.append({
            "Unsur": el,
            "Ar": ar,
            "Qu": qu,
            "n": n,
            "n×Ar": n_ar,
            "µc": muc,
            "µc²": muc2
        })

    u_bm = math.sqrt(sum_uc2)

    be = bm / n_equiv

    u_be = u_bm / n_equiv

    return bm, u_bm, be, u_be, rows, None

# ─────────────────────────────────────────────
# FUNGSI KETIDAKPASTIAN
# ─────────────────────────────────────────────
def mu_kalibrasi(u):
    return u / math.sqrt(3)

def mu_efek_temp(v, dt, alpha=2.1e-4):
    return v * (dt / math.sqrt(3)) * alpha

def mu_volume_gabungan(u_spek, volume, delta_T):

    mk = mu_kalibrasi(u_spek)

    met = mu_efek_temp(volume, delta_T)

    return math.sqrt(mk**2 + met**2)

# ─────────────────────────────────────────────
# SIDEBAR
# ─────────────────────────────────────────────
with st.sidebar:

    st.header("Pengaturan")

    T_spek = st.number_input(
        "Suhu spesifikasi alat (°C)",
        value=20.0
    )

    T_ruang = st.number_input(
        "Suhu ruang (°C)",
        value=31.0
    )

    delta_T = abs(T_ruang - T_spek)

    st.info(f"ΔT = {delta_T:.1f} °C")

# ─────────────────────────────────────────────
# TABS
# ─────────────────────────────────────────────
tab1, tab2 = st.tabs([
    "📥 Input",
    "📊 Hasil"
])

# ═══════════════════════════════════════════
# TAB INPUT
# ═══════════════════════════════════════════
with tab1:

    st.subheader("Input Data")

    formula = st.text_input(
        "Formula Senyawa",
        value="H2C2O4"
    )

    n_equiv = st.number_input(
        "Ekivalen per mol",
        value=2,
        min_value=1
    )

    u_buret = st.number_input(
        "U Spek Buret (mL)",
        value=0.05
    )

    vt1 = st.number_input(
        "VT Replikasi 1",
        value=27.50
    )

    vt2 = st.number_input(
        "VT Replikasi 2",
        value=27.53
    )

    Y = st.number_input(
        "Nilai Y",
        value=0.09093,
        format="%.6f"
    )

# ═══════════════════════════════════════════
# TAB HASIL
# ═══════════════════════════════════════════
with tab2:

    if st.button("🔬 Hitung", type="primary"):

        bm, u_bm, be, u_be, rows, err = calc_mu_BM(
            formula,
            n_equiv
        )

        if err:

            st.error(err)

        else:

            vt_rata = (vt1 + vt2) / 2

            mu_vt = mu_volume_gabungan(
                u_buret,
                vt_rata,
                delta_T
            )

            rel_be = u_be / be

            rel_vt = mu_vt / vt_rata

            uc_rel = math.sqrt(
                rel_be**2 + rel_vt**2
            )

            uc = Y * uc_rel

            U95 = 2 * uc

            c1, c2, c3, c4 = st.columns(4)

            c1.metric(
                "BM",
                f"{bm:.5f}"
            )

            c2.metric(
                "µ_BE",
                f"{u_be:.4e}"
            )

            c3.metric(
                "µ_VT",
                f"{mu_vt:.4e}"
            )

            c4.metric(
                "U95",
                f"{U95:.4e}"
            )

            st.success(
                f"Hasil akhir: {Y:.6f} ± {U95:.4e}"
            )

            st.divider()

            st.subheader("Breakdown Unsur")

            df = pd.DataFrame(rows)

            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True
            )

            st.divider()

            st.subheader("Rumus")

            st.latex(
                r"\mu_{BM} = \sqrt{\sum \mu_c^2}"
            )

            st.latex(
                r"u_c(Y)=Y\sqrt{\sum \left(\frac{\mu_i}{x_i}\right)^2}"
            )

            st.latex(
                r"U_{95}=2u_c(Y)"
            )
