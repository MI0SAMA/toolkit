import pytest
from click.testing import CliRunner
from toolkit.cli.main import main
from toolkit.cli.convert import convert
from toolkit.cli.analyze import analyze
from toolkit.cli.manipulate import manipulate
from conftest import data_path


@pytest.fixture
def runner():
    return CliRunner()


class TestCLIMain:
    def test_main_help(self, runner):
        result = runner.invoke(main, ['--help'])
        assert result.exit_code == 0
        assert 'convert' in result.output
        assert 'analyze' in result.output
        assert 'manipulate' in result.output

    def test_version(self, runner):
        result = runner.invoke(main, ['--version'])
        assert result.exit_code == 0
        assert '0.1.0' in result.output


class TestConvertTovasp:
    def test_vasp_to_vasp(self, runner, tmp_path):
        outfile = str(tmp_path / 'out.vasp')
        result = runner.invoke(main, [
            'convert', 'tovasp', data_path('water.vasp'), '-o', outfile
        ])
        assert result.exit_code == 0
        assert 'Written' in result.output

    def test_xyz_to_vasp(self, runner, tmp_path):
        outfile = str(tmp_path / 'out.vasp')
        lattice = '10,0,0,0,10,0,0,0,10'
        result = runner.invoke(main, [
            'convert', 'tovasp', data_path('water.xyz'),
            '-o', outfile, '--lattice', lattice
        ])
        assert result.exit_code == 0

    def test_xyz_without_lattice_fails(self, runner, tmp_path):
        outfile = str(tmp_path / 'out.vasp')
        result = runner.invoke(main, [
            'convert', 'tovasp', data_path('water.xyz'), '-o', outfile
        ])
        assert result.exit_code != 0

    def test_cartesian_flag(self, runner, tmp_path):
        outfile = str(tmp_path / 'out.vasp')
        result = runner.invoke(main, [
            'convert', 'tovasp', data_path('water.vasp'),
            '-o', outfile, '--cartesian'
        ])
        assert result.exit_code == 0


class TestConvertC2D:
    def test_c2d_with_lattice(self, runner, tmp_path):
        outfile = str(tmp_path / 'out.vasp')
        lattice = '10,0,0,0,10,0,0,0,10'
        result = runner.invoke(main, [
            'convert', 'c2d', data_path('water.xyz'),
            '-o', outfile, '--lattice', lattice
        ])
        assert result.exit_code == 0


class TestAnalyzeDensity:
    def test_density_basic(self, runner, tmp_path):
        outfile = str(tmp_path / 'density.csv')
        lattice = '15,0,0,0,15,0,0,0,15'
        result = runner.invoke(main, [
            'analyze', 'density', data_path('water.xyz'),
            '--lattice', lattice, '--element', 'O',
            '--interval', '1.0', '--min', '0', '--max', '15',
            '-o', outfile, '--no-plot'
        ])
        assert result.exit_code == 0
        assert 'Processed 1 files' in result.output


class TestAnalyzeAngle:
    def test_angle_basic(self, runner, tmp_path):
        lattice = '30,0,0,0,30,0,0,0,30'
        result = runner.invoke(main, [
            'analyze', 'angle', data_path('water.xyz'),
            '--lattice', lattice, '--period', 'xy',
            '--main', 'O', '--sub', 'H', '--cutoff', '1.5',
            '-o', str(tmp_path / 'angle')
        ])
        # May or may not find bonds depending on structure, but shouldn't crash
        assert result.exit_code == 0


class TestManipulateMove:
    def test_move_atoms(self, runner, tmp_path):
        outfile = str(tmp_path / 'moved.vasp')
        result = runner.invoke(main, [
            'manipulate', 'move', data_path('water.vasp'),
            '--atoms', '1', '--displacement', '0.5 0 0',
            '-o', outfile
        ])
        assert result.exit_code == 0

    def test_move_missing_args_prompts(self, runner, tmp_path):
        outfile = str(tmp_path / 'moved.vasp')
        result = runner.invoke(main, [
            'manipulate', 'move', data_path('water.vasp'),
            '-o', outfile
        ], input='1\n0.5 0 0\n')
        assert result.exit_code == 0


class TestManipulateJoint:
    def test_joint_basic(self, runner, tmp_path):
        outfile = str(tmp_path / 'joint.vasp')
        result = runner.invoke(main, [
            'manipulate', 'joint',
            data_path('water.vasp'), data_path('water.vasp'),
            '--axis', 'z', '--gap', '0.5',
            '-o', outfile
        ])
        assert result.exit_code == 0


class TestManipulateExtend:
    def test_extend_basic(self, runner, tmp_path):
        outfile = str(tmp_path / 'ext.vasp')
        result = runner.invoke(main, [
            'manipulate', 'extend', data_path('water.vasp'),
            '--nx', '2', '--ny', '2', '--nz', '1',
            '-o', outfile
        ])
        assert result.exit_code == 0
