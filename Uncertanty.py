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

    # support dot notation e.g. Na2B4O7.10H2O
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
        rows.append({"Unsur": el, "Ar": ar, "Qu": qu, "n": n, "n×Ar": n_ar, "µc": muc, "µc²": muc2})
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
    mk = mu_kalibrasi(u_spek)
    met = mu_efek_temp(volume, delta_T, alpha)
    mep = mu_end_point(v_nst, delta_T, alpha)
    return math.sqrt(mk**2 + met**2 + mep**2), mk, met, mep

def mu_neraca(u_neraca, k=math.sqrt(3)):
    return (2 * u_neraca) / k

def mu_kemurnian(u_purity, k=math.sqrt(3)):
    return u_purity / k

def mu_faktor_pengali(v_lt, u_lt, v_pip, u_pip, delta_T, fp, alpha=2.1e-4):
    mu_lt = math.sqrt(mu_kalibrasi(u_lt)**2 + mu_efek_temp(v_lt, delta_T, alpha)**2)
    mu_pip = math.sqrt(mu_kalibrasi(u_pip)**2 + mu_efek_temp(v_pip, delta_T, alpha)**2)
    return fp * math.sqrt((mu_lt/v_lt)**2 + (mu_pip/v_pip)**2), mu_lt, mu_pip

def mu_presisi(std_dev, n_rep):
    return std_dev / math.sqrt(n_rep)

def mu_normalitas_from_U95(U95):
    return U95 / 2.0


# ─────────────────────────────────────────────
# UI STREAMLIT
# ─────────────────────────────────────────────
st.set_page_config(
    page_title="Estimasi Ketidakpastian Titrimetri",
    page_icon="⚗️",
    layout="wide"
)

st.title("⚗️ Estimasi Ketidakpastian Titrimetri")
st.caption("Kalkulasi µ komponen dan ketidakpastian gabungan U₉₅ secara otomatis")

# ── SIDEBAR: Pilih jenis titrasi ──────────────────────────────────
with st.sidebar:
    st.header("Pengaturan")
    jenis = st.selectbox("Jenis Analisis", [
        "Standarisasi larutan (1 baku primer)",
        "Penetapan kadar analit",
    ])

    st.divider()
    st.subheader("Pilih Komponen µ")
    use_neraca    = st.checkbox("µ Bobot baku primer (neraca)", value=True)
    use_be        = st.checkbox("µ Bobot Ekivalen (µ_BE)", value=True)
    use_vt        = st.checkbox("µ Volume Titran (buret)", value=True)
    use_fp        = st.checkbox("µ Faktor Pengali (labu + pipet)", value=True)
    use_purity    = st.checkbox("µ Kemurnian baku primer", value=True)
    use_vs        = st.checkbox("µ Volume Sampel (pipet sampel)", value=False)
    use_n_std     = st.checkbox("µ Normalitas larutan standar (dari U₉₅)", value=False)
    use_pm        = st.checkbox("µ Presisi Metode (repeatability)", value=False)

    st.divider()
    st.subheader("Data Sekunder Alat")
    alpha = st.number_input("Koef. muai air α (°C⁻¹)", value=2.1e-4, format="%.2e")
    T_spek = st.number_input("Suhu spesifikasi alat (°C)", value=20.0)
    T_ruang = st.number_input("Suhu ruang (°C)", value=31.0)
    delta_T = abs(T_ruang - T_spek)
    st.info(f"ΔT = {delta_T:.1f} °C")

# ── MAIN CONTENT ──────────────────────────────────────────────────
tabs = st.tabs(["📥 Input Data", "📊 Hasil & Breakdown", "📐 Kalkulator µ_BM/BE"])

# ═══════════════════════════════════════════
# TAB 1 — INPUT DATA
# ═══════════════════════════════════════════
with tabs[0]:
    inputs = {}

    col1, col2 = st.columns(2)

    # ── BOBOT BAKU PRIMER ──────────────────
    if use_neraca:
        with col1:
            st.subheader("⚖️ Bobot Baku Primer")
            inputs["u_neraca"] = st.number_input("U neraca (mg)", value=0.2, format="%.4f",
                help="Ketidakpastian diperluas dari sertifikat kalibrasi neraca")
            inputs["bobot_primer"] = st.number_input("Bobot baku primer (mg)", value=630.5, format="%.4f")

    # ── BOBOT EKIVALEN ─────────────────────
    if use_be:
        with col2:
            st.subheader("🧪 Bobot Ekivalen Senyawa")
            inputs["formula_primer"] = st.text_input("Formula senyawa baku primer",
                value="H2C2O4", help="Contoh: NaOH, H2C2O4, Na2B4O7")
            inputs["n_equiv_primer"] = st.number_input("Ekivalen per mol (n)", value=2, min_value=1, step=1,
                help="n untuk BE = BM/n. Contoh: NaOH=1, H2C2O4=2, K2Cr2O7=6")

    # ── VOLUME TITRAN ─────────────────────
    if use_vt:
        st.subheader("🧫 Volume Titran (Buret)")
        c1, c2, c3, c4 = st.columns(4)
        inputs["u_buret"]   = c1.number_input("U spek buret (mL)", value=0.05, format="%.4f")
        inputs["vt_1"]      = c2.number_input("Volume titran rep-1 (mL)", value=27.50, format="%.4f")
        inputs["vt_2"]      = c3.number_input("Volume titran rep-2 (mL)", value=27.53, format="%.4f")
        inputs["v_nst"]     = c4.number_input("½ skala terkecil buret (mL)", value=0.005, format="%.4f",
            help="Biasanya 0.005 mL untuk buret 50 mL")
        vt_rata = (inputs["vt_1"] + inputs["vt_2"]) / 2
        st.info(f"Rata-rata volume titran: **{vt_rata:.4f} mL**")
        inputs["vt_rata"] = vt_rata

    # ── FAKTOR PENGALI ────────────────────
    if use_fp:
        st.subheader("🔢 Faktor Pengali (Labu Takar + Pipet)")
        c1, c2, c3, c4, c5 = st.columns(5)
        inputs["v_lt"]  = c1.number_input("Volume labu takar (mL)", value=100.0, format="%.2f")
        inputs["u_lt"]  = c2.number_input("U spek labu takar (mL)", value=0.2, format="%.4f")
        inputs["v_pip"] = c3.number_input("Volume pipet (mL)", value=25.0, format="%.2f")
        inputs["u_pip"] = c4.number_input("U spek pipet (mL)", value=0.03, format="%.4f")
        inputs["fp"]    = c5.number_input("Faktor pengali FP (VLT/Vpip)", value=4.0, format="%.4f")

    # ── KEMURNIAN ─────────────────────────
    if use_purity:
        st.subheader("💊 Kemurnian Baku Primer")
        c1, c2 = st.columns(2)
        inputs["purity"]   = c1.number_input("Kemurnian (%)", value=99.0, format="%.2f")
        inputs["u_purity"] = c2.number_input("U kemurnian (%)", value=1.0, format="%.4f")

    # ── VOLUME SAMPEL ─────────────────────
    if use_vs:
        st.subheader("🧴 Volume Sampel")
        c1, c2 = st.columns(2)
        inputs["v_sampel"] = c1.number_input("Volume sampel (mL)", value=25.0, format="%.2f")
        inputs["u_vs_pip"] = c2.number_input("U spek pipet sampel (mL)", value=0.03, format="%.4f")

    # ── NORMALITAS LARUTAN STANDAR ────────
    if use_n_std:
        st.subheader("📌 Normalitas Larutan Standar")
        c1, c2 = st.columns(2)
        inputs["N_std"]    = c1.number_input("Normalitas larutan standar (mgrek/mL)", value=0.09093, format="%.6f")
        inputs["U95_nstd"] = c2.number_input("U₉₅ normalitas standar", value=2.25e-4, format="%.2e",
            help="Ambil dari hasil standarisasi sebelumnya")

    # ── PRESISI METODE ────────────────────
    if use_pm:
        st.subheader("📈 Presisi Metode (Repeatability)")
        c1, c2, c3 = st.columns(3)
        kadar_1 = c1.number_input("Kadar hasil rep-1 (%)", value=0.303, format="%.6f")
        kadar_2 = c2.number_input("Kadar hasil rep-2 (%)", value=0.305, format="%.6f")
        inputs["n_rep"] = c3.number_input("Jumlah replikasi (n)", value=2, min_value=2, step=1)
        kadar_list = [kadar_1, kadar_2]
        mean_k = sum(kadar_list) / len(kadar_list)
        std_k  = math.sqrt(sum((x - mean_k)**2 for x in kadar_list) / (len(kadar_list) - 1))
        inputs["std_kadar"] = std_k
        inputs["mean_kadar"] = mean_k
        st.info(f"Rerata kadar: {mean_k:.6f}% | Std Deviasi: {std_k:.6f}%")

    # ── NILAI Y (HASIL PENETAPAN) ─────────
    st.divider()
    st.subheader("📌 Nilai Hasil Pengukuran (Y)")
    inputs["Y"] = st.number_input(
        "Nilai Y (normalitas / kadar yang ditetapkan)",
        value=0.09093,
        format="%.6f",
        help="Masukkan nilai normalitas atau kadar yang dihasilkan dari titrasi"
    )
    inputs["satuan_Y"] = st.text_input("Satuan Y", value="mgrek/mL")


# ═══════════════════════════════════════════
# TAB 2 — HASIL & BREAKDOWN
# ═══════════════════════════════════════════
with tabs[1]:
    st.subheader("Hasil Perhitungan Ketidakpastian")

    if st.button("🔬 Hitung Ketidakpastian", type="primary"):

        components = []   # list of dict {simbol, uraian, xi, satuan, mu_xi, mu_rel, mu_rel2}
        errors = []

        # 1. µ NERACA
        if use_neraca:
            um = mu_neraca(inputs["u_neraca"])
            m  = inputs["bobot_primer"]
            components.append({
                "Simbol": "µm",
                "Uraian": "Bobot baku primer (neraca)",
                "Xi": m,
                "Satuan": "mg",
                "µ(Xi)": um,
                "µ(Xi)/Xi": um / m if m else 0,
                "(µ(Xi)/Xi)²": (um / m)**2 if m else 0,
                "Detail": f"µm = 2×U_neraca/√3 = 2×{inputs['u_neraca']}/√3 = {um:.6f} mg"
            })

        # 2. µ BE
        if use_be:
            formula = inputs.get("formula_primer", "")
            n_eq    = inputs.get("n_equiv_primer", 1)
            bm, u_bm, be, u_be, rows_be, err = calc_mu_BM(formula, n_eq)
            if err:
                errors.append(f"µ_BE: {err}")
            else:
                components.append({
                    "Simbol": "µBE",
                    "Uraian": f"Bobot Ekivalen {formula}",
                    "Xi": be,
                    "Satuan": "mg/mgrek",
                    "µ(Xi)": u_be,
                    "µ(Xi)/Xi": u_be / be if be else 0,
                    "(µ(Xi)/Xi)²": (u_be / be)**2 if be else 0,
                    "Detail": f"BM={bm:.5f} g/mol | µ_BM={u_bm:.4e} | BE=BM/{n_eq}={be:.5f} | µ_BE={u_be:.4e}"
                })

        # 3. µ VOLUME TITRAN
        if use_vt:
            vt   = inputs["vt_rata"]
            mu_vt, mk, met, mep = mu_volume_gabungan(
                inputs["u_buret"], vt, delta_T, inputs["v_nst"], alpha
            )
            components.append({
                "Simbol": "µVT",
                "Uraian": "Volume titran (buret)",
                "Xi": vt,
                "Satuan": "mL",
                "µ(Xi)": mu_vt,
                "µ(Xi)/Xi": mu_vt / vt if vt else 0,
                "(µ(Xi)/Xi)²": (mu_vt / vt)**2 if vt else 0,
                "Detail": f"µ_kal={mk:.4e} | µ_ET={met:.4e} | µ_EP={mep:.4e} | µ_VT=√Σ={mu_vt:.4e}"
            })

        # 4. µ FAKTOR PENGALI
        if use_fp:
            mu_fp, mu_lt, mu_pip = mu_faktor_pengali(
                inputs["v_lt"], inputs["u_lt"],
                inputs["v_pip"], inputs["u_pip"],
                delta_T, inputs["fp"], alpha
            )
            fp = inputs["fp"]
            components.append({
                "Simbol": "µFP",
                "Uraian": "Faktor pengali (labu + pipet)",
                "Xi": fp,
                "Satuan": "",
                "µ(Xi)": mu_fp,
                "µ(Xi)/Xi": mu_fp / fp if fp else 0,
                "(µ(Xi)/Xi)²": (mu_fp / fp)**2 if fp else 0,
                "Detail": f"µ_LT={mu_lt:.4e} mL | µ_pip={mu_pip:.4e} mL | µ_FP={mu_fp:.4e}"
            })

        # 5. µ KEMURNIAN
        if use_purity:
            mp = mu_kemurnian(inputs["u_purity"])
            p  = inputs["purity"]
            components.append({
                "Simbol": "µP",
                "Uraian": "Kemurnian baku primer",
                "Xi": p,
                "Satuan": "%",
                "µ(Xi)": mp,
                "µ(Xi)/Xi": mp / p if p else 0,
                "(µ(Xi)/Xi)²": (mp / p)**2 if p else 0,
                "Detail": f"µP = U_purity/√3 = {inputs['u_purity']}/√3 = {mp:.6f} %"
            })

        # 6. µ VOLUME SAMPEL
        if use_vs:
            vs    = inputs["v_sampel"]
            mu_vs = math.sqrt(
                mu_kalibrasi(inputs["u_vs_pip"])**2 +
                mu_efek_temp(vs, delta_T, alpha)**2
            )
            components.append({
                "Simbol": "µVS",
                "Uraian": "Volume sampel (pipet)",
                "Xi": vs,
                "Satuan": "mL",
                "µ(Xi)": mu_vs,
                "µ(Xi)/Xi": mu_vs / vs if vs else 0,
                "(µ(Xi)/Xi)²": (mu_vs / vs)**2 if vs else 0,
                "Detail": f"µVS = √(µ_kal² + µ_ET²) = {mu_vs:.4e} mL"
            })

        # 7. µ NORMALITAS STANDAR
        if use_n_std:
            mu_n = mu_normalitas_from_U95(inputs["U95_nstd"])
            N    = inputs["N_std"]
            components.append({
                "Simbol": "µN",
                "Uraian": "Normalitas larutan standar",
                "Xi": N,
                "Satuan": "mgrek/mL",
                "µ(Xi)": mu_n,
                "µ(Xi)/Xi": mu_n / N if N else 0,
                "(µ(Xi)/Xi)²": (mu_n / N)**2 if N else 0,
                "Detail": f"µN = U₉₅/2 = {inputs['U95_nstd']:.4e}/2 = {mu_n:.4e}"
            })

        # 8. µ PRESISI METODE
        if use_pm:
            mu_pm = mu_presisi(inputs["std_kadar"], inputs["n_rep"])
            mean_k = inputs["mean_kadar"]
            components.append({
                "Simbol": "µPM",
                "Uraian": "Presisi metode (repeatability)",
                "Xi": mean_k,
                "Satuan": "%",
                "µ(Xi)": mu_pm,
                "µ(Xi)/Xi": mu_pm / mean_k if mean_k else 0,
                "(µ(Xi)/Xi)²": (mu_pm / mean_k)**2 if mean_k else 0,
                "Detail": f"µPM = s/√n = {inputs['std_kadar']:.6f}/√{inputs['n_rep']} = {mu_pm:.6f}"
            })

        # ── TAMPILKAN ERROR ──────────────────
        for e in errors:
            st.error(e)

        if components:
            Y = inputs["Y"]
            sat_Y = inputs["satuan_Y"]

            sum_rel2 = sum(c["(µ(Xi)/Xi)²"] for c in components)
            sqrt_sum = math.sqrt(sum_rel2)
            u_c      = Y * sqrt_sum
            U95      = 2 * u_c

            # ── METRIC CARDS ─────────────────
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Y (nilai penetapan)", f"{Y:.6g}", sat_Y)
            m2.metric("√Σ(µᵢ/xᵢ)²", f"{sqrt_sum:.6e}")
            m3.metric("u_c(Y)", f"{u_c:.6e}", sat_Y)
            m4.metric("U₉₅ (k=2)", f"{U95:.6e}", sat_Y)

            st.success(f"**Hasil: Y = {Y:.6g} ± {U95:.4e} {sat_Y}**")

            st.divider()

            # ── TABEL BREAKDOWN ──────────────
            st.subheader("Tabel Breakdown Komponen Ketidakpastian")
            df = pd.DataFrame(components)[
                ["Simbol","Uraian","Xi","Satuan","µ(Xi)","µ(Xi)/Xi","(µ(Xi)/Xi)²"]
            ]

            # Format scientific notation
            for col in ["µ(Xi)","µ(Xi)/Xi","(µ(Xi)/Xi)²"]:
                df[col] = df[col].apply(lambda x: f"{x:.4e}")

            st.dataframe(df, use_container_width=True, hide_index=True)

            # ── RUMUS DETAIL PER KOMPONEN ────
            st.divider()
            st.subheader("Detail Perhitungan per Komponen")
            for c in components:
                with st.expander(f"**{c['Simbol']}** — {c['Uraian']}"):
                    st.code(c["Detail"])
                    st.markdown(f"""
| Parameter | Nilai |
|-----------|-------|
| Xᵢ | {c['Xi']:.6g} {c['Satuan']} |
| µ(Xᵢ) | {c['µ(Xi)']:.4e} |
| µ(Xᵢ)/Xᵢ | {c['µ(Xi)/Xi']:.4e} |
| (µ(Xᵢ)/Xᵢ)² | {c['(µ(Xi)/Xi)²']:.4e} |
""")

            # ── RUMUS PROPAGASI ──────────────
            st.divider()
            st.subheader("Rumus Propagasi Gabungan")
            st.latex(r"u_c(Y) = Y \cdot \sqrt{\sum_i \left(\frac{\mu_i}{x_i}\right)^2}")
            st.latex(r"U_{95} = 2 \times u_c(Y)")
            st.markdown(f"""
| Langkah | Nilai |
|---------|-------|
| Σ(µᵢ/xᵢ)² | `{sum_rel2:.6e}` |
| √Σ(µᵢ/xᵢ)² | `{sqrt_sum:.6e}` |
| Y | `{Y:.6g}` {sat_Y} |
| u_c(Y) = Y × √Σ | `{u_c:.6e}` |
| **U₉₅ = 2 × u_c** | **`{U95:.6e}`** |
""")

        else:
            st.warning("Pilih minimal satu komponen µ di sidebar, lalu klik Hitung.")


# ═══════════════════════════════════════════
# TAB 3 — KALKULATOR µ_BM/BE
# ═══════════════════════════════════════════
with tabs[2]:
    st.subheader("🔬 Kalkulator µ_BM dan µ_BE")
    st.caption("Masukkan formula kimia → sistem hitung otomatis dari data tabel periodik IUPAC 2021")

    col1, col2 = st.columns([2, 1])
    with col1:
        formula_inp = st.text_input(
            "Formula senyawa",
            value="H2C2O4",
            placeholder="Contoh: NaOH, K2Cr2O7, Na2B4O7.10H2O"
        )
    with col2:
        n_equiv_inp = st.number_input(
            "Ekivalen per mol (n)",
            value=2, min_value=1, step=1,
            help="BE = BM / n"
        )

    st.caption("Contoh: NaOH · KMnO4 · Na2CO3 · H2C2O4 · Na2B4O7 · K2Cr2O7 · FeCl3 · CH3COOH")

    if st.button("Hitung µ_BM & µ_BE", type="primary"):
        bm, u_bm, be, u_be, rows_be, err = calc_mu_BM(formula_inp, n_equiv_inp)
        if err:
            st.error(err)
        else:
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("BM", f"{bm:.5f}", "g/mol")
            c2.metric("µ_BM", f"{u_bm:.4e}", "mg/mmol")
            c3.metric("BE = BM/n", f"{be:.5f}", "mg/mgrek")
            c4.metric("µ_BE", f"{u_be:.4e}", "mg/mgrek")

            st.divider()
            st.subheader("Breakdown per Unsur")
            df_be = pd.DataFrame(rows_be)
            for col in ["µc", "µc²"]:
                df_be[col] = df_be[col].apply(lambda x: f"{x:.4e}")
            df_be["n×Ar"] = df_be["n×Ar"].apply(lambda x: f"{x:.5f}")
            df_be["Ar"]   = df_be["Ar"].apply(lambda x: f"{x:.5f}")
            df_be["Qu"]   = df_be["Qu"].apply(lambda x: f"{x:.2e}")
            st.dataframe(df_be, use_container_width=True, hide_index=True)

            st.divider()
            st.subheader("Rumus")
            st.latex(r"\mu_{c,i} = \frac{n_i \times Q_{u,i}}{\sqrt{3}}")
            st.latex(r"\mu_{BM} = \sqrt{\sum_i \mu_{c,i}^2}")
            st.latex(r"\mu_{BE} = \frac{\mu_{BM}}{n_{ekuivalen}}")

    st.divider()
    st.subheader("📋 Database Unsur (IUPAC 2021)")
    db_df = pd.DataFrame([
        {"Unsur": el, "Ar (g/mol)": v["ar"], "Qu": v["qu"]}
        for el, v in PERIODIC_DB.items()
    ])
    db_df["Qu"] = db_df["Qu"].apply(lambda x: f"{x:.2e}")
    st.dataframe(db_df, use_container_width=True, hide_index=True, height=300)
