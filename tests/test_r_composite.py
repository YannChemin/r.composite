"""Tests for r.composite.

Default mode writes three true-colour FCELL maps (<output>.r/g/b).
Legacy palette mode (-p) writes a single indexed CELL map.
"""

import os
import pytest
import grass.script as gs


@pytest.fixture
def session(tmp_path):
    """GRASS session with a 3×3 test region and known RGB input maps."""
    project = tmp_path / "r_composite_project"
    gs.create_project(project)
    with gs.setup.init(project, env=os.environ.copy()) as sess:
        gs.run_command("g.region", n=3, s=0, e=3, w=0, res=1, env=sess.env)

        # red varies by row (50, 100, 150 N→S), green by column (50, 100, 150 W→E)
        gs.mapcalc("red   = row() * 50", overwrite=True, env=sess.env)
        gs.mapcalc("green = col() * 50", overwrite=True, env=sess.env)
        gs.mapcalc("blue  = 100",        overwrite=True, env=sess.env)

        # Apply explicit greyscale tables so colour-lookup values are predictable:
        # GRASS "grey" maps [min, max] → [black, white].
        # With range [50, 150]: value 50→R=0, 100→R=127, 150→R=255.
        for name in ("red", "green", "blue"):
            gs.run_command("r.colors", map=name, color="grey", env=sess.env)

        yield sess


# ─── True-colour (default) mode ──────────────────────────────────────────────

class TestTrueColor:
    @pytest.fixture(autouse=True)
    def run_composite(self, session):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="tc",
            env=session.env,
        )
        self.session = session

    def _ascii(self, name):
        raw = gs.read_command(
            "r.out.ascii", input=name, flags="h", env=self.session.env
        )
        return [[float(v) for v in row.split()] for row in raw.strip().splitlines()]

    def test_output_maps_are_fcell(self):
        for suffix in (".r", ".g", ".b"):
            info = gs.raster_info("tc" + suffix, env=self.session.env)
            assert info["datatype"] == "FCELL", \
                f"tc{suffix} should be FCELL, got {info['datatype']}"

    def test_red_constant_within_row(self):
        """Red component derives from row() so must be identical across columns."""
        rows = self._ascii("tc.r")
        for row in rows:
            assert len(set(row)) == 1, \
                f"Red values should be constant per row, got {row}"

    def test_red_increases_south(self):
        """Red component increases from north (row 1) to south (row 3)."""
        rows = self._ascii("tc.r")
        row_vals = [r[0] for r in rows]
        assert row_vals[0] < row_vals[1] < row_vals[2], \
            f"Red should increase N→S, got {row_vals}"

    def test_green_constant_within_column(self):
        """Green component derives from col() so must be identical across rows."""
        rows = self._ascii("tc.g")
        ncols = len(rows[0])
        for col_idx in range(ncols):
            col_vals = [rows[r][col_idx] for r in range(len(rows))]
            assert len(set(col_vals)) == 1, \
                f"Green values should be constant per column, got {col_vals}"

    def test_green_increases_east(self):
        """Green component increases from west (col 1) to east (col 3)."""
        rows = self._ascii("tc.g")
        col_vals = rows[0]
        assert col_vals[0] < col_vals[1] < col_vals[2], \
            f"Green should increase W→E, got {col_vals}"

    def test_blue_constant(self):
        """Blue component is constant (blue = 100)."""
        rows = self._ascii("tc.b")
        all_vals = [v for row in rows for v in row]
        assert len(set(all_vals)) == 1, \
            f"Blue should be constant, got {set(all_vals)}"

    def test_values_in_byte_range(self):
        """All component values must be in [0, 255]."""
        for suffix in (".r", ".g", ".b"):
            rows = self._ascii("tc" + suffix)
            for row in rows:
                for v in row:
                    assert 0.0 <= v <= 255.0, \
                        f"tc{suffix} value {v} out of [0, 255]"


# ─── NULL handling (true-colour) ─────────────────────────────────────────────

class TestNullHandling:
    @pytest.fixture(autouse=True)
    def setup(self, session):
        self.session = session
        gs.mapcalc(
            "red_null = if(row() == 2, null(), row() * 50)",
            overwrite=True, env=session.env,
        )
        gs.run_command("r.colors", map="red_null", color="grey", env=session.env)

    def test_null_propagation_to_red_band(self):
        """NULLs in the red input appear as NULLs in output.r."""
        gs.run_command(
            "r.composite",
            red="red_null", green="green", blue="blue",
            output="null_tc",
            env=self.session.env,
        )
        raw = gs.read_command(
            "r.out.ascii", input="null_tc.r", env=self.session.env
        )
        assert "*" in raw, "Expected NULL ('*') in null_tc.r"

    def test_null_isolated_to_red_band(self):
        """NULLs in the red input must NOT affect output.g or output.b."""
        gs.run_command(
            "r.composite",
            red="red_null", green="green", blue="blue",
            output="null_tc",
            env=self.session.env,
        )
        for suffix in (".g", ".b"):
            raw = gs.read_command(
                "r.out.ascii", input="null_tc" + suffix,
                env=self.session.env
            )
            assert "*" not in raw, \
                f"null_tc{suffix} should have no NULLs (only red band is null)"

    def test_null_fill_eliminates_nulls(self):
        """null_value= replaces NULLs with the fill value instead of propagating."""
        gs.run_command(
            "r.composite",
            red="red_null", green="green", blue="blue",
            output="filled_tc",
            null_value="128",
            env=self.session.env,
        )
        for suffix in (".r", ".g", ".b"):
            raw = gs.read_command(
                "r.out.ascii", input="filled_tc" + suffix,
                env=self.session.env
            )
            assert "*" not in raw, \
                f"filled_tc{suffix} should have no NULLs with null_value=128"


# ─── Palette (legacy -p) mode ────────────────────────────────────────────────

class TestPaletteMode:
    @pytest.fixture(autouse=True)
    def run_palette(self, session):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="pal",
            flags="p",
            env=session.env,
        )
        self.session = session

    def _cell_rows(self, name):
        raw = gs.read_command(
            "r.out.ascii", input=name, flags="h", env=self.session.env
        )
        return [[int(v) for v in row.split()] for row in raw.strip().splitlines()]

    def test_output_is_cell(self):
        info = gs.raster_info("pal", env=self.session.env)
        assert info["datatype"] == "CELL"

    def test_single_output_map(self):
        """Palette mode produces exactly one output map (not .r/.g/.b)."""
        result = gs.find_file("pal", element="raster", env=self.session.env)
        assert result["name"] == "pal"
        for suffix in (".r", ".g", ".b"):
            missing = gs.find_file("pal" + suffix, element="raster",
                                   env=self.session.env)
            assert missing["name"] == "", \
                f"pal{suffix} should not exist in palette mode"

    def test_grid_size(self):
        rows = self._cell_rows("pal")
        assert len(rows) == 3
        for row in rows:
            assert len(row) == 3

    def test_red_varies_by_row(self):
        """Palette index must change between rows (red channel varies with row)."""
        rows = self._cell_rows("pal")
        # All columns in the same row should be identical (red and blue are
        # constant per row; green changes per column but that changes the index)
        # so just verify that row 1 ≠ row 3
        assert rows[0][0] != rows[2][0], \
            "Palette index should differ between row 1 and row 3"

    def test_palette_null_propagation(self):
        """NULLs in any input propagate to the palette output."""
        gs.mapcalc(
            "red_null_pal = if(row() == 2, null(), row() * 50)",
            overwrite=True, env=self.session.env,
        )
        gs.run_command("r.colors", map="red_null_pal", color="grey",
                       env=self.session.env)
        gs.run_command(
            "r.composite",
            red="red_null_pal", green="green", blue="blue",
            output="pal_null",
            flags="p",
            env=self.session.env,
        )
        raw = gs.read_command(
            "r.out.ascii", input="pal_null", env=self.session.env
        )
        assert "*" in raw, "Expected NULL ('*') in palette output"


# ─── Palette dithering modes ──────────────────────────────────────────────────

class TestPaletteDithering:
    """Floyd-Steinberg and Bayer dithering produce valid CELL output."""

    @pytest.fixture(autouse=True)
    def setup(self, session):
        self.session = session

    def _cell_rows(self, name):
        raw = gs.read_command(
            "r.out.ascii", input=name, flags="h", env=self.session.env
        )
        return [[int(v) for v in row.split()] for row in raw.strip().splitlines()]

    def test_floyd_steinberg(self):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="pal_fs",
            flags="pd",
            env=self.session.env,
        )
        rows = self._cell_rows("pal_fs")
        assert len(rows) == 3 and all(len(r) == 3 for r in rows)

    def test_bayer(self):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="pal_bayer",
            flags="pb",
            env=self.session.env,
        )
        rows = self._cell_rows("pal_bayer")
        assert len(rows) == 3 and all(len(r) == 3 for r in rows)

    def test_gamma(self):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="pal_gamma",
            flags="pg",
            env=self.session.env,
        )
        info = gs.raster_info("pal_gamma", env=self.session.env)
        assert info["datatype"] == "CELL"

    def test_oklab(self):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="pal_oklab",
            flags="pk",
            env=self.session.env,
        )
        info = gs.raster_info("pal_oklab", env=self.session.env)
        assert info["datatype"] == "CELL"


# ─── Adaptive levels ─────────────────────────────────────────────────────────

class TestAdaptiveLevels:
    def test_auto_levels_produces_valid_output(self, session):
        gs.run_command(
            "r.composite",
            red="red", green="green", blue="blue",
            output="pal_auto",
            levels="auto",
            flags="p",
            env=session.env,
        )
        info = gs.raster_info("pal_auto", env=session.env)
        assert info["datatype"] == "CELL"
