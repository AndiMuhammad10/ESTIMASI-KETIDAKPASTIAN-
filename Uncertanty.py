import streamlit as st
import math
import pandas as pd

# ─────────────────────────────────────────────
# DATABASE TABEL PERIODIK (Ar + Qu IUPAC 2021)
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
    "Cr": {"ar": 51.9961,   "qu": 6e-4},
    "Ba": {"ar": 137.327,   "qu": 7e-4},
    "Sn": {"ar": 118.710,   "qu": 7e-3},
    "Pb": {"ar": 207.2,     "qu": 1.1},
    "Hg": {"ar": 200.59,    "qu": 4e-2},
    "Ti": {"ar": 47.867,    "qu": 1e-3},
    "V":  {"ar": 50.9415,   "qu": 1e-5},
    "As": {"ar": 74.9216,   "qu": 2e-5},
    "Se": {"ar": 78.96,     "qu": 3e-3},
}

# ─────────────────────────────────────────────
# PARSER FORMULA KIMIA
# ─────────────────────────────────────────────
def parse_formula(formula: str) -> dict:
    import re
    result = {}

    def parse_group(s, mult=1):
        counts = {}
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)')
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


def calc_mu_BM(formula: str, n_equiv: int = 1):
    parsed = parse_formula(formula)
    sqrt3 = math.sqrt(3)
    bm = 0.0
    sum_uc2 = 0.0
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
        rows.append({"Unsur": el, "Ar": ar, "Qu": qu, "n": n,
                     "n×Ar": n_ar, "µc": muc, "µc²": muc2})
    u_bm = math.sqrt(sum_uc2)
    be = bm / n_equiv
    u_be = u_bm / n_equiv
    return bm, u_bm, be, u_be, rows, None


# ─────────────────────────────────────────────
# FUNGSI HITUNG KOMPONEN µ
# ─────────────────────────────────────────────
def mu_kalibrasi(u_spek):
    return u_spek / math.sqrt(3)

def mu_efek_temp(volume, delta_T, alpha=2.1e-4):
    return volume * (delta_T / math.sqrt(3)) * alpha

def mu_end_point(v_nst, delta_T, alpha=2.1e-4):
    return v_nst * (delta_T / math.sqrt(3)) * alpha

def mu_volume_gabungan(u_spek, volume, delta_T, v_nst=0.005, alpha=2.1e-4):
    mk  = mu_kalibrasi(u_spek)
    met = mu_efek_temp(volume, delta_T, alpha)
    mep = mu_end_point(v_nst, delta_T, alpha)
    return math.sqrt(mk**2 + met**2 + mep**2), mk, met, mep

def mu_neraca(u_neraca):
    return (2 * u_neraca) / math.sqrt(3)

def mu_kemurnian(u_purity):
    return u_purity / math.sqrt(3)

def mu_faktor_pengali(v_lt, u_lt, v_pip, u_pip, delta_T, fp, alpha=2.1e-4):
    mu_lt  = math.sqrt(mu_kalibrasi(u_lt)**2  + mu_efek_temp(v_lt,  delta_T, alpha)**2)
    mu_pip = math.sqrt(mu_kalibrasi(u_pip)**2 + mu_efek_temp(v_pip, delta_T, alpha)**2)
    return fp * math.sqrt((mu_lt/v_lt)**2 + (mu_pip/v_pip)**2), mu_lt, mu_pip

def mu_presisi(std_dev, n_rep):
    return std_dev / math.sqrt(n_rep)

def mu_normalitas_from_U95(U95):
    return U95 / 2.0


# ─────────────────────────────────────────────
# HELPER: KOTAK PERHITUNGAN LANGKAH-DEMI-LANGKAH
# ─────────────────────────────────────────────
def show_step_box(label: str, steps: list[str], result_line: str):
    """
    Tampilkan expander berisi langkah perhitungan per komponen µ.
    steps  : list string tiap langkah
    result : baris hasil akhir
    """
    with st.expander(f"📐 Cara perhitungan {label}", expanded=False):
        for i, s in enumerate(steps, 1):
            st.markdown(f"**Langkah {i}:** {s}")
        st.success(f"**Hasil → {result_line}**")


# ═══════════════════════════════════════════
# STREAMLIT CONFIG
# ═══════════════════════════════════════════
st.set_page_config(
    page_title="Estimasi Ketidakpastian Titrimetri",
    page_icon="⚗️",
    layout="wide"
)

st.title("⚗️ Estimasi Ketidakpastian Titrimetri")
st.caption("Kalkulasi µ komponen dan ketidakpastian gabungan U₉₅ secara otomatis")

# ── SIDEBAR ────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Pengaturan")
    jenis = st.selectbox("Jenis Analisis", [
        "Standarisasi larutan (1 baku primer)",
        "Penetapan kadar analit",
    ])
    st.divider()
    st.subheader("Pilih Komponen µ")
    use_neraca = st.checkbox("µ Bobot baku primer (neraca)", value=True)
    use_be     = st.checkbox("µ Bobot Ekivalen (µ_BE)",     value=True)
    use_vt     = st.checkbox("µ Volume Titran (buret)",      value=True)
    use_fp     = st.checkbox("µ Faktor Pengali (labu+pipet)",value=True)
    use_purity = st.checkbox("µ Kemurnian baku primer",      value=True)
    use_vs     = st.checkbox("µ Volume Sampel (pipet sampel)",value=False)
    use_n_std  = st.checkbox("µ Normalitas standar (U₉₅)",  value=False)
    use_pm     = st.checkbox("µ Presisi Metode (repeatability)", value=False)

    st.divider()
    st.subheader("Data Sekunder Alat")
    alpha   = st.number_input("Koef. muai air α (°C⁻¹)", value=2.1e-4, format="%.2e")
    T_spek  = st.number_input("Suhu spesifikasi alat (°C)", value=20.0)
    T_ruang = st.number_input("Suhu ruang (°C)", value=31.0)
    delta_T = abs(T_ruang - T_spek)
    st.info(f"ΔT = {delta_T:.1f} °C")


# ── TABS ───────────────────────────────────────────────────────────
tabs = st.tabs(["📥 Input & Perhitungan µ", "📊 Hasil Akhir U₉₅", "📐 Kalkulator µ_BM/BE"])

inputs = {}

# ═══════════════════════════════════════════
# TAB 1 — INPUT + TAMPILAN PERHITUNGAN µ
# ═══════════════════════════════════════════
with tabs[0]:

    # ──────────────────────────────────────
    # 1. µ NERACA
    # ──────────────────────────────────────
    if use_neraca:
        st.markdown("---")
        st.subheader("⚖️ 1. Ketidakpastian Bobot Baku Primer (µm)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["u_neraca"]      = st.number_input("U neraca sertifikat (mg)", value=0.2, format="%.4f",
                help="Ketidakpastian diperluas dari sertifikat kalibrasi neraca")
            inputs["bobot_primer"]  = st.number_input("Bobot baku primer ditimbang (mg)", value=630.5, format="%.4f")

        um = mu_neraca(inputs["u_neraca"])
        inputs["_mu_neraca"] = um

        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_m = \frac{2 \times U_{neraca}}{\sqrt{3}}")
        show_step_box("µm (bobot baku primer)", [
            f"Penimbangan dilakukan **2 kali** (tara + isi), sehingga U total = 2 × U_neraca",
            f"Distribusi **rectangular** → k = √3",
            f"µm = 2 × {inputs['u_neraca']:.4f} / √3",
            f"µm = {2*inputs['u_neraca']:.4f} / {math.sqrt(3):.6f}",
        ], f"µm = **{um:.6f} mg**")


    # ──────────────────────────────────────
    # 2. µ BOBOT EKIVALEN
    # ──────────────────────────────────────
    if use_be:
        st.markdown("---")
        st.subheader("🧪 2. Ketidakpastian Bobot Ekivalen (µ_BE)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["formula_primer"]  = st.text_input("Formula senyawa baku primer",
                value="H2C2O4", help="Contoh: NaOH, H2C2O4, Na2B4O7.10H2O")
            inputs["n_equiv_primer"]  = st.number_input("Ekivalen per mol (n)", value=2, min_value=1, step=1,
                help="BE = BM/n. NaOH=1, H₂C₂O₄=2, K₂Cr₂O₇=6")
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_{c,i} = \frac{n_i \times Q_{u,i}}{\sqrt{3}}")
            st.latex(r"\mu_{BM} = \sqrt{\sum_i \mu_{c,i}^2} \quad ; \quad \mu_{BE} = \frac{\mu_{BM}}{n}")

        formula = inputs["formula_primer"]
        n_eq    = inputs["n_equiv_primer"]
        bm, u_bm, be, u_be, rows_be, err_be = calc_mu_BM(formula, n_eq)

        if err_be:
            st.error(err_be)
            inputs["_mu_be"] = 0; inputs["_be_val"] = 1
        else:
            inputs["_mu_be"]  = u_be
            inputs["_be_val"] = be

            parsed = parse_formula(formula)
            sqrt3  = math.sqrt(3)

            steps_be = [f"Parse formula **{formula}** → unsur: " +
                        ", ".join([f"{el} (n={n})" for el, n in parsed.items()])]
            for el, n in parsed.items():
                ar = PERIODIC_DB[el]["ar"]
                qu = PERIODIC_DB[el]["qu"]
                muc = n * qu / sqrt3
                steps_be.append(
                    f"**{el}**: µc = {n} × {qu:.2e} / √3 = **{muc:.4e}** mg/mmol"
                )
            sum_uc2 = sum((parsed[el] * PERIODIC_DB[el]["qu"] / sqrt3)**2 for el in parsed)
            steps_be.append(f"Σµc² = {sum_uc2:.4e}")
            steps_be.append(f"µ_BM = √{sum_uc2:.4e} = **{u_bm:.4e}** mg/mmol")
            steps_be.append(f"BE = BM / n = {bm:.5f} / {n_eq} = **{be:.5f}** mg/mgrek")
            steps_be.append(f"µ_BE = µ_BM / n = {u_bm:.4e} / {n_eq} = **{u_be:.4e}** mg/mgrek")

            show_step_box("µ_BE (bobot ekivalen)", steps_be,
                          f"µ_BE = {u_be:.4e} mg/mgrek  |  BE = {be:.5f} mg/mgrek")

            df_be = pd.DataFrame(rows_be)
            for c in ["µc","µc²"]:
                df_be[c] = df_be[c].apply(lambda x: f"{x:.4e}")
            df_be["n×Ar"] = df_be["n×Ar"].apply(lambda x: f"{x:.5f}")
            df_be["Qu"]   = df_be["Qu"].apply(lambda x: f"{x:.2e}")
            st.dataframe(df_be, use_container_width=True, hide_index=True)


    # ──────────────────────────────────────
    # 3. µ VOLUME TITRAN
    # ──────────────────────────────────────
    if use_vt:
        st.markdown("---")
        st.subheader("🧫 3. Ketidakpastian Volume Titran (µ_VT)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["u_buret"] = st.number_input("U spek buret (mL)", value=0.05, format="%.4f")
            inputs["vt_1"]    = st.number_input("Volume titran rep-1 (mL)", value=27.50, format="%.4f")
            inputs["vt_2"]    = st.number_input("Volume titran rep-2 (mL)", value=27.53, format="%.4f")
            inputs["v_nst"]   = st.number_input("½ skala terkecil buret (mL)", value=0.005, format="%.4f")
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_{kal} = \frac{U_{spek}}{\sqrt{3}}")
            st.latex(r"\mu_{ET} = \frac{V \cdot \Delta T}{\sqrt{3}} \cdot \alpha")
            st.latex(r"\mu_{EP} = \frac{V_{NST} \cdot \Delta T}{\sqrt{3}} \cdot \alpha")
            st.latex(r"\mu_{VT} = \sqrt{\mu_{kal}^2 + \mu_{ET}^2 + \mu_{EP}^2}")

        vt_rata = (inputs["vt_1"] + inputs["vt_2"]) / 2
        inputs["vt_rata"] = vt_rata
        st.info(f"Rata-rata VT = ({inputs['vt_1']} + {inputs['vt_2']}) / 2 = **{vt_rata:.4f} mL**")

        mu_vt, mk, met, mep = mu_volume_gabungan(
            inputs["u_buret"], vt_rata, delta_T, inputs["v_nst"], alpha)
        inputs["_mu_vt"] = mu_vt

        show_step_box("µ_VT (volume titran)", [
            f"**Sub-komponen 1 — Kalibrasi buret (Tipe B)**",
            f"µ_kal = U_spek / √3 = {inputs['u_buret']:.4f} / √3 = **{mk:.6f} mL**",
            f"**Sub-komponen 2 — Efek Temperatur (Tipe B)**",
            f"ΔT = |{T_ruang:.1f} − {T_spek:.1f}| = {delta_T:.1f} °C | α = {alpha:.2e} °C⁻¹",
            f"µ_ET = V × (ΔT/√3) × α = {vt_rata:.4f} × ({delta_T:.1f}/√3) × {alpha:.2e} = **{met:.6f} mL**",
            f"**Sub-komponen 3 — Titik Akhir / End Point (Tipe A)**",
            f"µ_EP = V_NST × (ΔT/√3) × α = {inputs['v_nst']:.4f} × ({delta_T:.1f}/√3) × {alpha:.2e} = **{mep:.2e} mL**",
            f"**Gabungan:**",
            f"µ_VT = √({mk:.4e}² + {met:.4e}² + {mep:.4e}²)",
        ], f"µ_VT = **{mu_vt:.6f} mL**")


    # ──────────────────────────────────────
    # 4. µ FAKTOR PENGALI
    # ──────────────────────────────────────
    if use_fp:
        st.markdown("---")
        st.subheader("🔢 4. Ketidakpastian Faktor Pengali (µ_FP)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["v_lt"]  = st.number_input("Volume labu takar (mL)", value=100.0, format="%.2f")
            inputs["u_lt"]  = st.number_input("U spek labu takar (mL)", value=0.2,   format="%.4f")
            inputs["v_pip"] = st.number_input("Volume pipet (mL)",       value=25.0,  format="%.2f")
            inputs["u_pip"] = st.number_input("U spek pipet (mL)",       value=0.03,  format="%.4f")
            inputs["fp"]    = st.number_input("Faktor pengali FP (VLT/Vpipet)", value=4.0, format="%.4f")
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_{LT} = \sqrt{\left(\frac{U_{LT}}{\sqrt{3}}\right)^2 + \left(\frac{V_{LT} \cdot \Delta T \cdot \alpha}{\sqrt{3}}\right)^2}")
            st.latex(r"\mu_{pip} = \sqrt{\left(\frac{U_{pip}}{\sqrt{3}}\right)^2 + \left(\frac{V_{pip} \cdot \Delta T \cdot \alpha}{\sqrt{3}}\right)^2}")
            st.latex(r"\mu_{FP} = FP \cdot \sqrt{\left(\frac{\mu_{LT}}{V_{LT}}\right)^2 + \left(\frac{\mu_{pip}}{V_{pip}}\right)^2}")

        mu_fp, mu_lt, mu_pip_val = mu_faktor_pengali(
            inputs["v_lt"], inputs["u_lt"],
            inputs["v_pip"], inputs["u_pip"],
            delta_T, inputs["fp"], alpha)
        inputs["_mu_fp"] = mu_fp

        mk_lt  = mu_kalibrasi(inputs["u_lt"])
        met_lt = mu_efek_temp(inputs["v_lt"], delta_T, alpha)
        mk_pp  = mu_kalibrasi(inputs["u_pip"])
        met_pp = mu_efek_temp(inputs["v_pip"], delta_T, alpha)

        show_step_box("µ_FP (faktor pengali)", [
            f"**Labu Takar {inputs['v_lt']:.0f} mL:**",
            f"µ_kal_LT = {inputs['u_lt']:.4f}/√3 = **{mk_lt:.4e} mL**",
            f"µ_ET_LT  = {inputs['v_lt']:.2f} × ({delta_T:.1f}/√3) × {alpha:.2e} = **{met_lt:.4e} mL**",
            f"µ_LT = √({mk_lt:.4e}² + {met_lt:.4e}²) = **{mu_lt:.4e} mL**",
            f"**Pipet {inputs['v_pip']:.0f} mL:**",
            f"µ_kal_pip = {inputs['u_pip']:.4f}/√3 = **{mk_pp:.4e} mL**",
            f"µ_ET_pip  = {inputs['v_pip']:.2f} × ({delta_T:.1f}/√3) × {alpha:.2e} = **{met_pp:.4e} mL**",
            f"µ_pip = √({mk_pp:.4e}² + {met_pp:.4e}²) = **{mu_pip_val:.4e} mL**",
            f"**Gabungan FP = {inputs['fp']:.4f}:**",
            f"µ_FP = {inputs['fp']:.4f} × √(({mu_lt:.4e}/{inputs['v_lt']:.2f})² + ({mu_pip_val:.4e}/{inputs['v_pip']:.2f})²)",
        ], f"µ_FP = **{mu_fp:.6f}**")


    # ──────────────────────────────────────
    # 5. µ KEMURNIAN
    # ──────────────────────────────────────
    if use_purity:
        st.markdown("---")
        st.subheader("💊 5. Ketidakpastian Kemurnian (µ_P)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["purity"]   = st.number_input("Kemurnian baku primer (%)", value=99.0, format="%.2f")
            inputs["u_purity"] = st.number_input("U kemurnian (%)", value=1.0, format="%.4f")
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_P = \frac{U_P}{\sqrt{3}}")

        mp = mu_kemurnian(inputs["u_purity"])
        inputs["_mu_purity"] = mp

        show_step_box("µ_P (kemurnian)", [
            f"Distribusi **rectangular** → k = √3",
            f"µ_P = U_P / √3 = {inputs['u_purity']:.4f} / √3 = {inputs['u_purity']:.4f} / {math.sqrt(3):.6f}",
        ], f"µ_P = **{mp:.6f} %**")


    # ──────────────────────────────────────
    # 6. µ VOLUME SAMPEL
    # ──────────────────────────────────────
    if use_vs:
        st.markdown("---")
        st.subheader("🧴 6. Ketidakpastian Volume Sampel (µ_VS)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["v_sampel"] = st.number_input("Volume sampel (mL)", value=25.0, format="%.2f")
            inputs["u_vs_pip"] = st.number_input("U spek pipet sampel (mL)", value=0.03, format="%.4f")
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_{VS} = \sqrt{\left(\frac{U_{pip}}{\sqrt{3}}\right)^2 + \left(\frac{V_S \cdot \Delta T \cdot \alpha}{\sqrt{3}}\right)^2}")

        mk_vs  = mu_kalibrasi(inputs["u_vs_pip"])
        met_vs = mu_efek_temp(inputs["v_sampel"], delta_T, alpha)
        mu_vs  = math.sqrt(mk_vs**2 + met_vs**2)
        inputs["_mu_vs"] = mu_vs

        show_step_box("µ_VS (volume sampel)", [
            f"**Sub-komponen 1 — Kalibrasi pipet sampel:**",
            f"µ_kal = {inputs['u_vs_pip']:.4f} / √3 = **{mk_vs:.4e} mL**",
            f"**Sub-komponen 2 — Efek Temperatur:**",
            f"µ_ET = {inputs['v_sampel']:.2f} × ({delta_T:.1f}/√3) × {alpha:.2e} = **{met_vs:.4e} mL**",
            f"**Gabungan:**",
            f"µ_VS = √({mk_vs:.4e}² + {met_vs:.4e}²)",
        ], f"µ_VS = **{mu_vs:.6f} mL**")


    # ──────────────────────────────────────
    # 7. µ NORMALITAS STANDAR
    # ──────────────────────────────────────
    if use_n_std:
        st.markdown("---")
        st.subheader("📌 7. Ketidakpastian Normalitas Larutan Standar (µ_N)")
        col1, col2 = st.columns(2)
        with col1:
            inputs["N_std"]    = st.number_input("Normalitas larutan standar (mgrek/mL)",
                value=0.09093, format="%.6f")
            inputs["U95_nstd"] = st.number_input("U₉₅ normalitas standar",
                value=2.25e-4, format="%.2e",
                help="Dari laporan standarisasi sebelumnya")
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"\mu_N = \frac{U_{95}}{k} = \frac{U_{95}}{2}")

        mu_n = mu_normalitas_from_U95(inputs["U95_nstd"])
        inputs["_mu_n_std"] = mu_n

        show_step_box("µ_N (normalitas standar)", [
            f"U₉₅ diambil dari hasil standarisasi = **{inputs['U95_nstd']:.4e}**",
            f"Coverage factor k = 2 (tingkat kepercayaan 95%)",
            f"µ_N = U₉₅ / k = {inputs['U95_nstd']:.4e} / 2",
        ], f"µ_N = **{mu_n:.4e} mgrek/mL**")


    # ──────────────────────────────────────
    # 8. µ PRESISI METODE
    # ──────────────────────────────────────
    if use_pm:
        st.markdown("---")
        st.subheader("📈 8. Ketidakpastian Presisi Metode (µ_PM)")
        col1, col2 = st.columns(2)
        with col1:
            kadar_1 = st.number_input("Kadar rep-1 (%)", value=0.303, format="%.6f")
            kadar_2 = st.number_input("Kadar rep-2 (%)", value=0.305, format="%.6f")
            inputs["n_rep"] = st.number_input("Jumlah replikasi (n)", value=2, min_value=2, step=1)
        with col2:
            st.markdown("**Rumus:**")
            st.latex(r"s = \sqrt{\frac{\sum(x_i - \bar{x})^2}{n-1}}")
            st.latex(r"\mu_{PM} = \frac{s}{\sqrt{n}}")

        kadar_list  = [kadar_1, kadar_2]
        mean_k = sum(kadar_list) / len(kadar_list)
        std_k  = math.sqrt(sum((x - mean_k)**2 for x in kadar_list) / (len(kadar_list) - 1))
        mu_pm  = mu_presisi(std_k, inputs["n_rep"])
        inputs["std_kadar"]  = std_k
        inputs["mean_kadar"] = mean_k
        inputs["_mu_pm"]     = mu_pm

        dev1 = kadar_1 - mean_k
        dev2 = kadar_2 - mean_k
        show_step_box("µ_PM (presisi metode)", [
            f"Rerata x̄ = ({kadar_1:.6f} + {kadar_2:.6f}) / 2 = **{mean_k:.6f} %**",
            f"Deviasi: (x₁ - x̄) = {dev1:.6f} | (x₂ - x̄) = {dev2:.6f}",
            f"Σ(xᵢ-x̄)² = {dev1**2:.4e} + {dev2**2:.4e} = {dev1**2+dev2**2:.4e}",
            f"s = √(Σ(xᵢ-x̄)² / (n-1)) = √({dev1**2+dev2**2:.4e} / {len(kadar_list)-1}) = **{std_k:.6f} %**",
            f"µ_PM = s / √n = {std_k:.6f} / √{inputs['n_rep']} = {std_k:.6f} / {math.sqrt(inputs['n_rep']):.4f}",
        ], f"µ_PM = **{mu_pm:.6f} %**")

    # ── NILAI Y ──────────────────────────
    st.markdown("---")
    st.subheader("📌 Nilai Hasil Pengukuran (Y)")
    col1, col2 = st.columns(2)
    inputs["Y"]        = col1.number_input("Nilai Y (normalitas / kadar)", value=0.09093, format="%.6f")
    inputs["satuan_Y"] = col2.text_input("Satuan Y", value="mgrek/mL")


# ═══════════════════════════════════════════
# TAB 2 — HASIL AKHIR
# ═══════════════════════════════════════════
with tabs[1]:
    st.subheader("Hasil Ketidakpastian Gabungan U₉₅")

    if st.button("🔬 Hitung Ketidakpastian Gabungan", type="primary"):
        components = []
        errors = []

        if use_neraca:
            um = inputs.get("_mu_neraca", mu_neraca(inputs.get("u_neraca", 0.2)))
            m  = inputs.get("bobot_primer", 1)
            components.append({
                "Simbol": "µm", "Uraian": "Bobot baku primer",
                "Xi": m, "Satuan": "mg", "µ(Xi)": um,
                "µ(Xi)/Xi": um/m if m else 0,
                "(µ(Xi)/Xi)²": (um/m)**2 if m else 0,
            })

        if use_be:
            u_be = inputs.get("_mu_be", 0)
            be   = inputs.get("_be_val", 1)
            components.append({
                "Simbol": "µBE", "Uraian": f"Bobot Ekivalen {inputs.get('formula_primer','')}",
                "Xi": be, "Satuan": "mg/mgrek", "µ(Xi)": u_be,
                "µ(Xi)/Xi": u_be/be if be else 0,
                "(µ(Xi)/Xi)²": (u_be/be)**2 if be else 0,
            })

        if use_vt:
            vt     = inputs.get("vt_rata", 1)
            mu_vt  = inputs.get("_mu_vt", 0)
            components.append({
                "Simbol": "µVT", "Uraian": "Volume titran (buret)",
                "Xi": vt, "Satuan": "mL", "µ(Xi)": mu_vt,
                "µ(Xi)/Xi": mu_vt/vt if vt else 0,
                "(µ(Xi)/Xi)²": (mu_vt/vt)**2 if vt else 0,
            })

        if use_fp:
            fp    = inputs.get("fp", 1)
            mu_fp = inputs.get("_mu_fp", 0)
            components.append({
                "Simbol": "µFP", "Uraian": "Faktor pengali",
                "Xi": fp, "Satuan": "", "µ(Xi)": mu_fp,
                "µ(Xi)/Xi": mu_fp/fp if fp else 0,
                "(µ(Xi)/Xi)²": (mu_fp/fp)**2 if fp else 0,
            })

        if use_purity:
            p  = inputs.get("purity", 1)
            mp = inputs.get("_mu_purity", mu_kemurnian(inputs.get("u_purity", 1)))
            components.append({
                "Simbol": "µP", "Uraian": "Kemurnian baku primer",
                "Xi": p, "Satuan": "%", "µ(Xi)": mp,
                "µ(Xi)/Xi": mp/p if p else 0,
                "(µ(Xi)/Xi)²": (mp/p)**2 if p else 0,
            })

        if use_vs:
            vs    = inputs.get("v_sampel", 1)
            mu_vs = inputs.get("_mu_vs", 0)
            components.append({
                "Simbol": "µVS", "Uraian": "Volume sampel",
                "Xi": vs, "Satuan": "mL", "µ(Xi)": mu_vs,
                "µ(Xi)/Xi": mu_vs/vs if vs else 0,
                "(µ(Xi)/Xi)²": (mu_vs/vs)**2 if vs else 0,
            })

        if use_n_std:
            N    = inputs.get("N_std", 1)
            mu_n = inputs.get("_mu_n_std", 0)
            components.append({
                "Simbol": "µN", "Uraian": "Normalitas larutan standar",
                "Xi": N, "Satuan": "mgrek/mL", "µ(Xi)": mu_n,
                "µ(Xi)/Xi": mu_n/N if N else 0,
                "(µ(Xi)/Xi)²": (mu_n/N)**2 if N else 0,
            })

        if use_pm:
            mk   = inputs.get("mean_kadar", 1)
            mu_pm= inputs.get("_mu_pm", 0)
            components.append({
                "Simbol": "µPM", "Uraian": "Presisi metode",
                "Xi": mk, "Satuan": "%", "µ(Xi)": mu_pm,
                "µ(Xi)/Xi": mu_pm/mk if mk else 0,
                "(µ(Xi)/Xi)²": (mu_pm/mk)**2 if mk else 0,
            })

        for e in errors:
            st.error(e)

        if components:
            Y      = inputs.get("Y", 1)
            sat_Y  = inputs.get("satuan_Y", "")
            sum_r2 = sum(c["(µ(Xi)/Xi)²"] for c in components)
            sqrt_s = math.sqrt(sum_r2)
            u_c    = Y * sqrt_s
            U95    = 2 * u_c

            # metric cards
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Y", f"{Y:.6g}", sat_Y)
            m2.metric("√Σ(µᵢ/xᵢ)²", f"{sqrt_s:.4e}")
            m3.metric("u_c(Y)", f"{u_c:.4e}", sat_Y)
            m4.metric("U₉₅ (k=2)", f"{U95:.4e}", sat_Y)

            st.success(f"### Hasil Akhir: Y = {Y:.6g} ± {U95:.4e} {sat_Y}")
            st.divider()

            # tabel breakdown
            st.subheader("Tabel Propagasi Ketidakpastian")
            df = pd.DataFrame(components)[
                ["Simbol","Uraian","Xi","Satuan","µ(Xi)","µ(Xi)/Xi","(µ(Xi)/Xi)²"]]
            for col in ["µ(Xi)","µ(Xi)/Xi","(µ(Xi)/Xi)²"]:
                df[col] = df[col].apply(lambda x: f"{x:.4e}")
            st.dataframe(df, use_container_width=True, hide_index=True)

            # baris sigma
            st.markdown(f"""
| | |
|---|---|
| **Σ(µᵢ/xᵢ)²** | `{sum_r2:.6e}` |
| **√Σ(µᵢ/xᵢ)²** | `{sqrt_s:.6e}` |
| **u_c(Y) = Y × √Σ** | `{u_c:.6e}` |
| **U₉₅ = 2 × u_c** | **`{U95:.6e}`** |
""")

            st.divider()
            st.subheader("Rumus Propagasi")
            st.latex(r"u_c(Y) = Y \cdot \sqrt{\sum_i \left(\frac{\mu_i}{x_i}\right)^2}")
            st.latex(r"U_{95} = 2 \times u_c(Y) \quad (k=2,\ \text{tingkat kepercayaan } 95\%)")

        else:
            st.warning("Pilih minimal satu komponen µ di sidebar.")


# ═══════════════════════════════════════════
# TAB 3 — KALKULATOR µ_BM/BE STANDALONE
# ═══════════════════════════════════════════
with tabs[2]:
    st.subheader("🔬 Kalkulator µ_BM dan µ_BE")
    st.caption("Masukkan formula kimia → hitung otomatis dari data tabel periodik IUPAC 2021")

    col1, col2 = st.columns([2, 1])
    with col1:
        formula_inp = st.text_input("Formula senyawa",
            value="H2C2O4", placeholder="NaOH, K2Cr2O7, Na2B4O7.10H2O",
            key="tab3_formula")
    with col2:
        n_equiv_inp = st.number_input("Ekivalen per mol (n)",
            value=2, min_value=1, step=1, key="tab3_nequiv")

    st.caption("Contoh: NaOH · KMnO4 · Na2CO3 · H2C2O4 · Na2B4O7 · K2Cr2O7 · FeCl3 · CH3COOH")

    if st.button("Hitung µ_BM & µ_BE", type="primary", key="tab3_btn"):
        bm, u_bm, be, u_be, rows_be, err = calc_mu_BM(formula_inp, n_equiv_inp)
        if err:
            st.error(err)
        else:
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("BM", f"{bm:.5f}", "g/mol")
            c2.metric("µ_BM", f"{u_bm:.4e}", "mg/mmol")
            c3.metric("BE = BM/n", f"{be:.5f}", "mg/mgrek")
            c4.metric("µ_BE", f"{u_be:.4e}", "mg/mgrek")

            st.divider()
            st.subheader("Breakdown per Unsur")
            df_be2 = pd.DataFrame(rows_be)
            for c in ["µc","µc²"]:
                df_be2[c] = df_be2[c].apply(lambda x: f"{x:.4e}")
            df_be2["n×Ar"] = df_be2["n×Ar"].apply(lambda x: f"{x:.5f}")
            df_be2["Qu"]   = df_be2["Qu"].apply(lambda x: f"{x:.2e}")
            st.dataframe(df_be2, use_container_width=True, hide_index=True)

            st.divider()
            parsed2 = parse_formula(formula_inp)
            sqrt3   = math.sqrt(3)
            st.subheader("Langkah Perhitungan Detail")
            for el, n in parsed2.items():
                ar  = PERIODIC_DB[el]["ar"]
                qu  = PERIODIC_DB[el]["qu"]
                muc = n * qu / sqrt3
                st.markdown(f"- **{el}** (n={n}): µc = {n} × {qu:.2e} / √3 = **{muc:.4e}** mg/mmol")
            sum_uc2_2 = sum((parsed2[el]*PERIODIC_DB[el]["qu"]/sqrt3)**2 for el in parsed2)
            st.markdown(f"- **Σµc²** = {sum_uc2_2:.4e}")
            st.markdown(f"- **µ_BM** = √{sum_uc2_2:.4e} = **{u_bm:.4e}** mg/mmol")
            st.markdown(f"- **µ_BE** = {u_bm:.4e} / {n_equiv_inp} = **{u_be:.4e}** mg/mgrek")

            st.subheader("Rumus")
            st.latex(r"\mu_{c,i} = \frac{n_i \times Q_{u,i}}{\sqrt{3}}")
            st.latex(r"\mu_{BM} = \sqrt{\sum_i \mu_{c,i}^2} \qquad \mu_{BE} = \frac{\mu_{BM}}{n}")

    st.divider()
    st.subheader("📋 Database Unsur (IUPAC 2021)")
    db_df = pd.DataFrame([
        {"Unsur": el, "Ar (g/mol)": v["ar"], "Qu": f"{v['qu']:.2e}"}
        for el, v in PERIODIC_DB.items()
    ])
    st.dataframe(db_df, use_container_width=True, hide_index=True, height=300)
