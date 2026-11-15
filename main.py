from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pymatgen.io.cif import CifParser
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import tempfile, os

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def health_check():
    return {"status": "ok", "engine": "pymatgen"}

@app.post("/calculate")
async def calc_xrd(file: UploadFile = File(...), wavelength: str = Form("CuKa")):
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as tmp:
        tmp.write(await file.read())
        tmp_path = tmp.name
    try:
        parser = CifParser(tmp_path, occupancy_tolerance=100.0)
        structure = parser.parse_structures(primitive=False)[0]

        try:
            sga = SpacegroupAnalyzer(structure, symprec=1e-3)
            sg_symbol = sga.get_space_group_symbol()
            sg_number = sga.get_space_group_number()
        except Exception:
            sg_symbol = "Unknown"
            sg_number = "-"

        calculator = XRDCalculator(wavelength=wavelength)
        pattern = calculator.get_pattern(structure, two_theta_range=(5, 100))

        peaks = [
            {
                "two_theta": round(float(x), 4),
                "intensity": round(float(y), 2),
                "hkl": [list(h["hkl"]) for h in hkls],
                "d": round(float(d), 5)
            }
            for x, y, hkls, d in zip(pattern.x, pattern.y, pattern.hkls, pattern.d_hkls)
        ]
        return {
            "status": "ok",
            "formula": structure.composition.reduced_formula,
            "space_group": f"{sg_symbol} (No.{sg_number})",
            "site_count": len(structure),
            "peaks": peaks
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
