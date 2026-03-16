from __future__ import annotations

import math
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lpmv


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "outputs"
REAL_MODEL_URL = (
    "https://icgem.gfz-potsdam.de/getseries/03_other/HUST/HUST-Grace2020/60x60/"
    "HUST-Grace2020-n60-200301.gfc"
)
REAL_MODEL_PATH = DATA_DIR / "HUST-Grace2020-n60-200301.gfc"


@dataclass
class HarmonicModel:
    c: np.ndarray
    s: np.ndarray
    normalization: str
    name: str

    @property
    def nmax(self) -> int:
        return self.c.shape[0] - 1


def fully_normalized_pnm(n: int, m: int, x: np.ndarray) -> np.ndarray:
    delta = 1 if m == 0 else 0
    factor = math.sqrt(
        (2 * n + 1) * (2 - delta) * math.factorial(n - m) / math.factorial(n + m)
    )
    return factor * lpmv(m, n, x)


def schmidt_normalized_pnm(n: int, m: int, x: np.ndarray) -> np.ndarray:
    delta = 1 if m == 0 else 0
    factor = math.sqrt((2 - delta) * math.factorial(n - m) / math.factorial(n + m))
    return factor * lpmv(m, n, x)


def normalized_pnm(n: int, m: int, x: np.ndarray, normalization: str) -> np.ndarray:
    if normalization == "fully":
        return fully_normalized_pnm(n, m, x)
    if normalization == "schmidt":
        return schmidt_normalized_pnm(n, m, x)
    raise ValueError(f"Normalizacao desconhecida: {normalization}")


def build_basis(
    theta: np.ndarray,
    lam: np.ndarray,
    n: int,
    m: int,
    trig: str,
    normalization: str = "fully",
) -> np.ndarray:
    x = np.cos(theta)
    pnm = normalized_pnm(n, m, x, normalization)
    if m == 0:
        return pnm
    if trig == "cos":
        return pnm * np.cos(m * lam)
    if trig == "sin":
        return pnm * np.sin(m * lam)
    raise ValueError(f"Componente trigonometrica desconhecida: {trig}")


def synthesize_field(
    theta: np.ndarray,
    lam: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    normalization: str = "fully",
    n_trunc: int | None = None,
) -> np.ndarray:
    nmax = c.shape[0] - 1 if n_trunc is None else min(n_trunc, c.shape[0] - 1)
    field = np.zeros_like(theta, dtype=float)
    for n in range(nmax + 1):
        for m in range(n + 1):
            field += c[n, m] * build_basis(theta, lam, n, m, "cos", normalization)
            if m > 0:
                field += s[n, m] * build_basis(theta, lam, n, m, "sin", normalization)
    return field


def integrate_coefficients_regular(
    theta: np.ndarray,
    lam: np.ndarray,
    field: np.ndarray,
    nmax: int,
    normalization: str = "fully",
) -> tuple[np.ndarray, np.ndarray]:
    c = np.zeros((nmax + 1, nmax + 1))
    s = np.zeros((nmax + 1, nmax + 1))

    weights = np.sin(theta)
    dtheta = np.diff(theta[:, 0]).mean()
    dlam = np.diff(lam[0, :]).mean()

    for n in range(nmax + 1):
        for m in range(n + 1):
            cos_basis = build_basis(theta, lam, n, m, "cos", normalization)
            cos_norm = np.sum(cos_basis * cos_basis * weights) * dtheta * dlam
            c[n, m] = np.sum(field * cos_basis * weights) * dtheta * dlam / cos_norm
            if m > 0:
                sin_basis = build_basis(theta, lam, n, m, "sin", normalization)
                sin_norm = np.sum(sin_basis * sin_basis * weights) * dtheta * dlam
                s[n, m] = (
                    np.sum(field * sin_basis * weights) * dtheta * dlam / sin_norm
                )
    return c, s


def enumerate_terms(nmax: int) -> list[tuple[int, int, str]]:
    terms: list[tuple[int, int, str]] = []
    for n in range(nmax + 1):
        for m in range(n + 1):
            terms.append((n, m, "cos"))
            if m > 0:
                terms.append((n, m, "sin"))
    return terms


def build_design_matrix(
    theta: np.ndarray,
    lam: np.ndarray,
    nmax: int,
    normalization: str = "fully",
) -> tuple[np.ndarray, list[tuple[int, int, str]]]:
    terms = enumerate_terms(nmax)
    matrix = np.zeros((theta.size, len(terms)))
    for column, (n, m, trig) in enumerate(terms):
        matrix[:, column] = build_basis(theta, lam, n, m, trig, normalization)
    return matrix, terms


def fit_spherical_harmonics_irregular(
    theta: np.ndarray,
    lam: np.ndarray,
    field: np.ndarray,
    nmax: int,
    normalization: str = "fully",
) -> tuple[np.ndarray, np.ndarray]:
    a, terms = build_design_matrix(theta, lam, nmax, normalization)
    solution, *_ = np.linalg.lstsq(a, field, rcond=None)
    c = np.zeros((nmax + 1, nmax + 1))
    s = np.zeros((nmax + 1, nmax + 1))
    for value, (n, m, trig) in zip(solution, terms):
        if trig == "cos":
            c[n, m] = value
        else:
            s[n, m] = value
    return c, s


def degree_power(c: np.ndarray, s: np.ndarray) -> np.ndarray:
    nmax = c.shape[0] - 1
    power = np.zeros(nmax + 1)
    for n in range(nmax + 1):
        for m in range(n + 1):
            power[n] += c[n, m] ** 2 + s[n, m] ** 2
    return power


def rmse(reference: np.ndarray, estimate: np.ndarray) -> float:
    return float(np.sqrt(np.mean((reference - estimate) ** 2)))


def save_global_map(
    theta: np.ndarray,
    lam: np.ndarray,
    values: np.ndarray,
    title: str,
    filename: str,
) -> None:
    lat = 90.0 - np.degrees(theta)
    lon = np.degrees(lam)
    fig, ax = plt.subplots(figsize=(11, 5))
    mesh = ax.pcolormesh(lon, lat, values, shading="auto", cmap="coolwarm")
    ax.set_title(title)
    ax.set_xlabel("Longitude (graus)")
    ax.set_ylabel("Latitude (graus)")
    fig.colorbar(mesh, ax=ax, label="Valor do campo")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=180)
    plt.close(fig)


def save_degree_power_plot(power: np.ndarray, title: str, filename: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    degrees = np.arange(power.size)
    ax.plot(degrees, power, marker="o", linewidth=1.5)
    ax.set_title(title)
    ax.set_xlabel("Grau n")
    ax.set_ylabel(r"$\sigma_n^2$")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=180)
    plt.close(fig)


def download_real_model() -> Path | None:
    DATA_DIR.mkdir(exist_ok=True)
    if REAL_MODEL_PATH.exists():
        return REAL_MODEL_PATH
    try:
        urllib.request.urlretrieve(REAL_MODEL_URL, REAL_MODEL_PATH)
        return REAL_MODEL_PATH
    except Exception:
        return None


def load_gfc_model(path: Path, nmax_limit: int = 20) -> HarmonicModel:
    rows: list[tuple[int, int, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 5 and parts[0].lower() == "gfc":
                n = int(parts[1])
                m = int(parts[2])
                if n <= nmax_limit:
                    rows.append((n, m, float(parts[3]), float(parts[4])))
    nmax = max(n for n, _, _, _ in rows)
    c = np.zeros((nmax + 1, nmax + 1))
    s = np.zeros((nmax + 1, nmax + 1))
    for n, m, cnm, snm in rows:
        c[n, m] = cnm
        s[n, m] = snm
    return HarmonicModel(c=c, s=s, normalization="fully", name=path.stem)


def format_coefficients(c: np.ndarray, s: np.ndarray, threshold: float = 1e-3) -> str:
    lines: list[str] = []
    nmax = c.shape[0] - 1
    for n in range(nmax + 1):
        for m in range(n + 1):
            if abs(c[n, m]) > threshold:
                lines.append(f"C[{n},{m}] = {c[n, m]: .6f}")
            if m > 0 and abs(s[n, m]) > threshold:
                lines.append(f"S[{n},{m}] = {s[n, m]: .6f}")
    return "\n".join(lines) if lines else "Nenhum coeficiente acima do limiar."


def write_report(lines: list[str]) -> None:
    (ROOT / "relatorio_resultados.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)

    report: list[str] = ["# Resultados do Projeto 1", ""]

    nmax = 6
    theta_vals = np.linspace(1e-3, math.pi - 1e-3, 91)
    lam_vals = np.linspace(0.0, 2 * math.pi, 181, endpoint=False)
    theta_grid, lam_grid = np.meshgrid(theta_vals, lam_vals, indexing="ij")

    c_true = np.zeros((nmax + 1, nmax + 1))
    s_true = np.zeros((nmax + 1, nmax + 1))
    c_true[2, 0] = 3.0
    c_true[3, 1] = 2.0
    s_true[4, 2] = -1.5

    field_true = synthesize_field(theta_grid, lam_grid, c_true, s_true)
    save_global_map(
        theta_grid,
        lam_grid,
        field_true,
        "Campo sintetico original",
        "mapa_campo_original.png",
    )

    c_regular, s_regular = integrate_coefficients_regular(
        theta_grid, lam_grid, field_true, nmax=4
    )
    field_regular = synthesize_field(theta_grid, lam_grid, c_regular, s_regular, n_trunc=4)
    save_global_map(
        theta_grid,
        lam_grid,
        field_regular,
        "Reconstrucao a partir da grade regular",
        "mapa_reconstrucao_regular.png",
    )
    save_degree_power_plot(
        degree_power(c_regular, s_regular),
        "Potencia espectral por grau - experimento 1",
        "potencia_grau_regular.png",
    )

    report.extend(
        [
            "## Objetivo 1: recuperar coeficientes por integracao numerica",
            "",
            "Coeficientes relevantes recuperados na grade regular:",
            "```text",
            format_coefficients(c_regular, s_regular),
            "```",
            f"RMSE da reconstrucao regular: {rmse(field_true, field_regular):.6e}",
            "",
        ]
    )

    rng = np.random.default_rng(42)
    point_count = 1500
    u = rng.uniform(-1.0, 1.0, point_count)
    theta_irregular = np.arccos(u)
    lam_irregular = rng.uniform(0.0, 2 * math.pi, point_count)
    field_irregular = synthesize_field(theta_irregular, lam_irregular, c_true, s_true)
    c_irregular, s_irregular = fit_spherical_harmonics_irregular(
        theta_irregular, lam_irregular, field_irregular, nmax=4
    )
    field_irregular_grid = synthesize_field(
        theta_grid, lam_grid, c_irregular, s_irregular, n_trunc=4
    )
    save_global_map(
        theta_grid,
        lam_grid,
        field_irregular_grid,
        "Reconstrucao a partir de pontos irregulares",
        "mapa_reconstrucao_irregular.png",
    )
    save_degree_power_plot(
        degree_power(c_irregular, s_irregular),
        "Potencia espectral por grau - experimento 2",
        "potencia_grau_irregular.png",
    )
    report.extend(
        [
            "## Objetivo 2: ajuste com observacoes irregulares",
            "",
            "Coeficientes relevantes ajustados por minimos quadrados:",
            "```text",
            format_coefficients(c_irregular, s_irregular),
            "```",
            (
                "RMSE da reconstrucao via pontos irregulares na grade global: "
                f"{rmse(field_true, field_irregular_grid):.6e}"
            ),
            "",
        ]
    )

    power_regular = degree_power(c_regular, s_regular)
    report.extend(
        [
            "## Objetivo 3: distribuicao espectral de energia por grau",
            "",
            "```text",
            "\n".join(
                f"grau {degree}: {value:.6f}" for degree, value in enumerate(power_regular)
            ),
            "```",
            "",
        ]
    )

    report.extend(["## Tarefas adicionais", ""])

    report.extend(
        [
            "### 1. Reconstrucao do campo",
            "",
            (
                "A reconstrucao foi executada para os coeficientes obtidos na grade regular e "
                "nos pontos irregulares. Os mapas foram salvos em `outputs/`."
            ),
            "",
        ]
    )

    report.append("### 2. Mapas globais truncados em diferentes graus")
    report.append("")
    for degree in (2, 3, 4):
        truncated = synthesize_field(theta_grid, lam_grid, c_regular, s_regular, n_trunc=degree)
        save_global_map(
            theta_grid,
            lam_grid,
            truncated,
            f"Mapa truncado ate grau {degree}",
            f"mapa_truncado_grau_{degree}.png",
        )
        report.append(f"- Mapa truncado em grau {degree} salvo em `outputs/mapa_truncado_grau_{degree}.png`.")
    report.append("")

    theta_sample = np.linspace(1e-3, math.pi - 1e-3, 400)
    x_sample = np.cos(theta_sample)
    p_full = fully_normalized_pnm(4, 2, x_sample)
    p_schmidt = schmidt_normalized_pnm(4, 2, x_sample)
    scale_ratio = float(np.max(np.abs(p_full)) / np.max(np.abs(p_schmidt)))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(np.degrees(theta_sample), p_full, label="Fully normalized")
    ax.plot(np.degrees(theta_sample), p_schmidt, label="Schmidt", linestyle="--")
    ax.set_xlabel("Colatitude (graus)")
    ax.set_ylabel(r"$\bar{P}_{42}$")
    ax.set_title("Comparacao entre normalizacoes")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "comparacao_normalizacoes.png", dpi=180)
    plt.close(fig)
    report.extend(
        [
            "### 3. Comparacao entre Schmidt e fully normalized",
            "",
            (
                "Foi gerado um grafico comparando `P42` nas duas normalizacoes. "
                f"A razao entre os picos absolutos foi {scale_ratio:.4f}."
            ),
            "",
        ]
    )

    report.append("### 4. Aplicacao a coeficientes geopotenciais reais")
    report.append("")
    real_model_path = download_real_model()
    if real_model_path is None:
        report.append(
            "Nao foi possivel baixar automaticamente o modelo real nesta execucao. "
            "O script continua pronto para processar um arquivo `.gfc` colocado em `data/`."
        )
        report.append("")
    else:
        real_model = load_gfc_model(real_model_path, nmax_limit=20)
        power_real = degree_power(real_model.c, real_model.s)
        save_degree_power_plot(
            power_real,
            f"Potencia espectral por grau - {real_model.name}",
            "potencia_grau_modelo_real.png",
        )
        field_real = synthesize_field(
            theta_grid, lam_grid, real_model.c, real_model.s, n_trunc=12
        )
        save_global_map(
            theta_grid,
            lam_grid,
            field_real,
            f"Campo do modelo real truncado ate grau 12 - {real_model.name}",
            "mapa_modelo_real_grau_12.png",
        )
        report.append(
            f"Modelo real processado: `{real_model_path.name}` com graus ate {real_model.nmax}."
        )
        report.append("Os produtos foram salvos em `outputs/potencia_grau_modelo_real.png` e `outputs/mapa_modelo_real_grau_12.png`.")
        report.append("")

    write_report(report)
    print("Analise concluida.")
    print(f"Relatorio salvo em: {ROOT / 'relatorio_resultados.md'}")
    print(f"Figuras salvas em: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
