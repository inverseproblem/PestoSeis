import numpy as np
from scipy import interpolate
from scipy import signal
import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pestoseis.seismicwaves2d
import unittest

# +
scale = 2
# scale = 4
dh = 1.0 / scale

# Domain characteristics
domain_size = np.array([100, 100]) * scale
nx = domain_size[0]
nz = domain_size[1]

# Sources and receivers
sources = np.array([50, 50]) * scale
receivers = np.array([[20, 20], [75, 25], [70, 70], [30, 65],])

# Material parameters
vp_ref = 3000.0
vs_ref = 2000.0
rho_ref = 2500.0

# Gradient material parameters
vp_ref_range = np.array([1500.0, 3000.0])
vs_ref_range = np.array([1000.0, 2000.0])
rho_ref_range = np.array([1250.0, 2500.0])

# Simulation length
t_sim = 0.075
dt = 0.00009
# dt = 0.00005
nt = int(t_sim / dt)
t = np.arange(0.0, nt * dt, dt)

# Source-time function
t0 = 0.003
f0 = 250.0

# Initiallize input parameter dictionary
inpar = {}
inpar["ntimesteps"] = nt
inpar["nx"] = nx
inpar["nz"] = nz
inpar["dt"] = dt
inpar["dh"] = dh
inpar["savesnapshot"] = False
inpar["boundcond"] = "PML"


# +
def import_pestoseis_data(filename: str) -> dict:
    """
    Imports an HDF5 receiver dataset generated using PestoSeis.

    Args:
        filename (str): Filename of the receiver file; must end in ``.h5``

    Returns:
        data (dict): Imported receiver data as a dictionary
    """
    data = {}
    with h5py.File(filename, "r") as f:
        if len(f["seism"].shape) == 3:
            # Get magnitude of vx and vz compontents and correct the sign
            data["waveforms"] = np.linalg.norm(f["seism"][:, :, :], axis=2)
            # Assign the correct +/- sign
            # mask = f["seism"][:, :, 1] < 0
            # data["waveforms"][mask] *= -1
        else:
            data["waveforms"] = f["seism"][:, :]
        data["nt"] = data["waveforms"].shape[1]
        data["start_time"] = 0.0
        data["end_time"] = f["dt"][()] * data["nt"]
        data["dh"] = float(f["dh"][()])
        data["nx"] = int(f["nx"][()])
        data["nz"] = int(f["nz"][()])
        data["recpos"] = f["recpos"][:, :]
        data["srcij"] = f["srcij"][:]

    data["t"] = np.linspace(
        data["start_time"], data["end_time"], data["waveforms"].shape[1]
    )
    return data


def import_salvus_data(filename: str, interp: bool = False) -> dict:
    """
    Imports an HDF5 receiver dataset generated using Salvus.

    Args:
        filename (str): Filename of the receiver file; must end in ``.h5``
        interp (bool, optional): Whether to interpolate the values to match the dimensions of the PestoSeis data

    Returns:
        data (dict): Imported receiver data as a dictionary
    """
    data = {}
    with h5py.File(filename, "r") as f:
        data["waveforms"] = f["waveforms"][:, :]
        data["nt"] = data["waveforms"].shape[1]
        data["start_time"] = f["start_time_in_seconds"][()]
        data["end_time"] = f["end_time_in_seconds"][()]

    data["t"] = np.linspace(
        data["start_time"], data["end_time"], data["waveforms"].shape[1]
    )

    if interp:
        interp_fun = interpolate.interp1d(data["t"], data["waveforms"])
        data["waveforms"] = interp_fun(t)
        data["t"] = t

    return data


def _plot_compare_pestoseis_salvus(
    pestoseis_file: str, salvus_file: str, plot_title: str = "PestoSeis vs Salvus"
) -> None:
    """
    Plot a comparison between a set of receiver data output from both PestoSeis and Salvus.

    Args:
        pestoseis_file (str): Name of the PestoSeis receiver file; must end in ``.h5``
        salvus_file (str): Name of the Salvus receiver file; must end in ``.h5``
        plot_title (str, optional): Optional plot title
    """
    if pestoseis_file[-3:] != ".h5" or salvus_file[-3:] != ".h5":
        raise ValueError("Input files must be of file type ``.h5``")

    # Import data
    pesto_data = import_pestoseis_data(pestoseis_file)
    salvus_data = import_salvus_data(salvus_file)

    # Plot the results
    n_recs = 4
    # _, ax = plt.subplots(nrows=n_recs + 1, ncols=1, figsize=[16, 10])
    fig = plt.figure(figsize=[16, 10])
    gs = fig.add_gridspec(n_recs, 2)
    ax = []

    colors = ["r", "b", "c", "g"]
    for i in range(n_recs):
        ax.append(fig.add_subplot(gs[i, 0]))
        ax[i].plot(
            salvus_data["t"], salvus_data["waveforms"][i, :], color="k", label="Salvus",
        )
        ax[i].plot(
            pesto_data["t"],
            pesto_data["waveforms"][i, :],
            color=colors[i],
            linestyle="--",
            label="PestoSeis",
        )

        ax[i].legend()
        ax[i].set_ylabel("Amplitude")

    ax[0].set_title(plot_title)
    ax[i].set_xlabel("Time [s]")

    # Add a plot to show the source-receiver configuration
    ax.append(fig.add_subplot(gs[:, 1]))
    markersize = 10
    for i, color in enumerate(colors):
        ax[-1].plot(
            pesto_data["recpos"][i, 0],
            pesto_data["recpos"][i, 1],
            markersize=markersize,
            marker="v",
            linestyle="",
            color=color,
        )

    ax[-1].plot(
        pesto_data["srcij"][0] * pesto_data["dh"],
        pesto_data["srcij"][1] * pesto_data["dh"],
        "k*",
        markersize=markersize,
        label="Source",
    )

    ax[-1].set_xlim([0.0, pesto_data["nx"] * pesto_data["dh"]])
    ax[-1].set_ylim([0.0, pesto_data["nz"] * pesto_data["dh"]])
    ax[-1].set_xlabel("x-Position [m]")
    ax[-1].set_ylabel("z-Position [m]")
    ax[-1].set_title("Source-Receiver Setup")
    ax[-1].set_aspect("equal")
    ax[-1].invert_yaxis()

    source_labels = Line2D(
        [0],
        [0],
        linestyle="",
        marker="*",
        color="k",
        markersize=markersize,
        label="Source",
    )
    receiver_labels = Line2D(
        [0],
        [0],
        linestyle="",
        marker="v",
        markerfacecolor="none",
        color="k",
        markersize=markersize,
        label="Receivers",
    )
    plt.legend(handles=[source_labels, receiver_labels])

    plt.show()


def _acoustic_ref_data(free_surface: bool, homogeneous: bool = True) -> None:
    """
    Function used for generating PestoSeis reference data for the acoustic case.

    The new reference file is placed in the ``pestoseis/tests`` directory.  In
    order to avoid accidentally overwriting old reference solutions, the user
    must manually move these to the ``pestoseis/tests/reference_solutions``
    directory if they wish to update the previous reference solutions.
    """
    if homogeneous:
        model = np.ones(domain_size, dtype=float) * vp_ref
    else:
        model = np.zeros(domain_size, dtype=float)
        for i in range(nx):
            model[i, :] = np.linspace(vp_ref_range[0], vp_ref_range[1], nz)

    inpar["freesurface"] = free_surface

    if homogeneous:
        suffix_1 = "hom"
    else:
        suffix_1 = "grad"

    if free_surface:
        suffix_2 = "free_surface"
    else:
        suffix_2 = "PML"

    stf = pestoseis.seismicwaves2d.sourcetimefuncs.rickersource(t, t0, f0)

    pestoseis.seismicwaves2d.acousticwaveprop2D.solveacoustic2D(
        inpar,
        sources,
        model,
        stf,
        f0,
        receivers,
        outfileh5=f"pestoseis_acoustic_{suffix_1}_{suffix_2}_ref.h5",
    )
    _plot_compare_pestoseis_salvus(
        f"pestoseis_acoustic_{suffix_1}_{suffix_2}_ref.h5",
        f"reference_solutions/salvus_acoustic_{suffix_1}_{suffix_2}_ref.h5",
    )

    data = import_pestoseis_data(f"pestoseis_acoustic_{suffix_1}_{suffix_2}_ref.h5")
    ref = import_salvus_data(
        f"reference_solutions/salvus_acoustic_{suffix_1}_{suffix_2}_ref.h5", interp=True
    )

    diff = np.linalg.norm(data["waveforms"] - ref["waveforms"])
    print(
        f"L2 Norm between Salvus and the new PestoSeis reference solution is {diff*100}\n"
    )


def _elastic_ref_data(free_surface: bool, homogeneous: bool = True) -> None:
    """
    Function used for generating PestoSeis reference data for the elastic case.

    The new reference file is placed in the ``pestoseis/tests`` directory.  In
    order to avoid accidentally overwriting old reference solutions, the user
    must manually move these to the ``pestoseis/tests/reference_solutions``
    directory if they wish to update the previous reference solutions.
    """
    if homogeneous:
        vp_model = np.ones(domain_size, dtype=float) * vp_ref
        vs_model = np.ones(domain_size, dtype=float) * vs_ref
        rho_model = np.ones(domain_size, dtype=float) * rho_ref
    else:
        vp_model = np.zeros(domain_size, dtype=float)
        vs_model = np.zeros(domain_size, dtype=float)
        rho_model = np.zeros(domain_size, dtype=float)
        for i in range(nx):
            vp_model[i, :] = np.linspace(vp_ref_range[0], vp_ref_range[1], nz)
            vs_model[i, :] = np.linspace(vs_ref_range[0], vs_ref_range[1], nz)
            rho_model[i, :] = np.linspace(rho_ref_range[0], rho_ref_range[1], nz)

    lambda_model = rho_model * (vp_model ** 2 - 2.0 * vs_model ** 2)
    mu_model = rho_model * vs_model ** 2

    rockprops = {}
    rockprops["lambda"] = lambda_model
    rockprops["mu"] = mu_model
    rockprops["rho"] = rho_model

    inpar["freesurface"] = free_surface

    if homogeneous:
        suffix_1 = "hom"
    else:
        suffix_1 = "grad"

    if free_surface:
        suffix_2 = "free_surface"
    else:
        suffix_2 = "PML"

    inpar["sourcetype"] = "MomTensor"
    inpar["momtensor"] = dict(xx=1.0, xz=0.0, zz=1.0)
    inpar["seismogrkind"] = "velocity"

    stf = pestoseis.seismicwaves2d.sourcetimefuncs.ricker_1st_derivative_source(
        t, t0, f0
    )

    pestoseis.seismicwaves2d.elasticwaveprop2D.solveelastic2D(
        inpar,
        rockprops,
        sources,
        stf,
        f0,
        receivers,
        outfileh5=f"pestoseis_elastic_{suffix_1}_{suffix_2}_ref.h5",
    )
    _plot_compare_pestoseis_salvus(
        f"pestoseis_elastic_{suffix_1}_{suffix_2}_ref.h5",
        f"reference_solutions/salvus_elastic_{suffix_1}_{suffix_2}_ref.h5",
    )

    data = import_pestoseis_data(f"pestoseis_elastic_{suffix_1}_{suffix_2}_ref.h5")
    ref = import_salvus_data(
        f"reference_solutions/salvus_elastic_{suffix_1}_{suffix_2}_ref.h5", interp=True
    )

    #     diff = np.linalg.norm(data["waveforms"] - ref["waveforms"])
    diff = np.sqrt(
        np.sum((data["waveforms"] - ref["waveforms"]) ** 2)
        / np.sum(data["waveforms"] ** 2)
    )
    print(
        f"L2 Norm between Salvus and the new PestoSeis reference solution is {diff}\n"
    )


# +
def _test_acoustic_homogeneous(free_surface: bool, threshold: float = 0.02) -> bool:
    """
    Test function for checking if the acoustic wave solver produces the expected
    results.

    Solves a simple forward problem of a wave propagating through a homogeneous
    medium.  The result is compared to a reference solution within the
    `reference_solutions` directory.  If the difference (measured in an L2
    sense) is greater than the user-defined threshold, then the test will fail.

    Args:
        free_surface (bool): Whether to treat the top surface of the domain as
        a free surface boundary condition.
        threshold (float, optional): Acceptable error threshold above which the
        test will fail.

    Returns:
        (bool): Whether the test passed or failed.
    """
    model = np.ones(domain_size, dtype=float) * vp_ref

    inpar["freesurface"] = free_surface

    if free_surface:
        suffix = "free_surface"
    else:
        suffix = "PML"

    stf = pestoseis.seismicwaves2d.sourcetimefuncs.rickersource(t, t0, f0)

    data, _ = pestoseis.seismicwaves2d.acousticwaveprop2D.solveacoustic2D(
        inpar,
        sources,
        model,
        stf,
        f0,
        receivers,
        outfileh5=f"pestoseis_acoustic_hom_{suffix}.h5",
    )

    ref = import_pestoseis_data(
        f"reference_solutions/pestoseis_acoustic_hom_{suffix}_ref.h5"
    )["waveforms"]

    # Ensure sampling between the reference data and the simulated data are the same
    ref_resampled = np.zeros_like(data)
    for i, trace in enumerate(ref):
        interp_fun = interpolate.interp1d(np.linspace(0, 1, ref.shape[1]), trace)
        ref_resampled[i, :] = interp_fun(np.linspace(0, 1, data.shape[1]))

    diff = np.mean(np.linalg.norm(data - ref_resampled, axis=1))
    print(
        f"L2 Norm between simulated data and the PestoSeis reference solution is {diff}\n"
    )

    return diff < threshold


def _test_acoustic_gradient(threshold: float = 0.02) -> bool:
    """
    Test function for checking if the acoustic wave solver produces the expected
    results.

    Solves a simple forward problem of a wave propagating through a medium with
    linearly varying material properties with depth.  The result is compared to
    a reference solution within the `reference_solutions` directory.  If the
    difference (measured in an L2 sense) is greater than the user-defined
    threshold, then the test will fail.

    Args:
        threshold (float, optional): Acceptable error threshold above which the
        test will fail.

    Returns:
        (bool): Whether the test passed or failed.
    """
    model = np.zeros(domain_size, dtype=float)
    for i in range(nx):
        model[i, :] = np.linspace(vp_ref_range[0], vp_ref_range[1], nz)

    inpar["freesurface"] = False
    suffix = "PML"

    stf = pestoseis.seismicwaves2d.sourcetimefuncs.rickersource(t, t0, f0)

    data, _ = pestoseis.seismicwaves2d.acousticwaveprop2D.solveacoustic2D(
        inpar,
        sources,
        model,
        stf,
        f0,
        receivers,
        outfileh5=f"pestoseis_acoustic_grad_{suffix}.h5",
    )

    ref = import_pestoseis_data(
        f"reference_solutions/pestoseis_acoustic_grad_{suffix}_ref.h5"
    )["waveforms"]

    # Ensure sampling between the reference data and the simulated data are the same
    ref_resampled = np.zeros_like(data)
    for i, trace in enumerate(ref):
        interp_fun = interpolate.interp1d(np.linspace(0, 1, ref.shape[1]), trace)
        ref_resampled[i, :] = interp_fun(np.linspace(0, 1, data.shape[1]))

    diff = np.mean(np.linalg.norm(data - ref_resampled, axis=1))
    print(
        f"L2 Norm between simulated data and the PestoSeis reference solution is {diff}\n"
    )

    return diff < threshold


def _test_elastic_homogeneous(free_surface: bool, threshold: float = 5e-9) -> bool:
    """
    Test function for checking if the elastic wave solver produces the expected
    results.

    Solves a simple forward problem of a wave propagating through a homogeneous
    medium.  The result is compared to a reference solution within the
    `reference_solutions` directory.  If the difference (measured in an L2
    sense) is greater than the user-defined threshold, then the test will fail.

    Args:
        free_surface (bool): Whether to treat the top surface of the domain as
        a free surface boundary condition.
        threshold (float, optional): Acceptable error threshold above which the
        test will fail.

    Returns:
        (bool): Whether the test passed or failed.
    """
    vp_model = np.ones(domain_size, dtype=float) * vp_ref
    vs_model = np.ones(domain_size, dtype=float) * vs_ref
    rho_model = np.ones(domain_size, dtype=float) * rho_ref

    lambda_model = rho_model * (vp_model ** 2 - 2.0 * vs_model ** 2)
    mu_model = rho_model * vs_model ** 2

    rockprops = {}
    rockprops["lambda"] = lambda_model
    rockprops["mu"] = mu_model
    rockprops["rho"] = rho_model

    inpar["freesurface"] = free_surface
    if free_surface:
        suffix = "free_surface"
    else:
        suffix = "PML"

    inpar["sourcetype"] = "MomTensor"
    inpar["momtensor"] = dict(xx=1.0, xz=0.0, zz=1.0)
    inpar["seismogrkind"] = "velocity"

    stf = pestoseis.seismicwaves2d.sourcetimefuncs.ricker_1st_derivative_source(
        t, t0, f0
    )

    data, _ = pestoseis.seismicwaves2d.elasticwaveprop2D.solveelastic2D(
        inpar,
        rockprops,
        sources,
        stf,
        f0,
        receivers,
        outfileh5=f"pestoseis_elastic_hom_{suffix}.h5",
    )

    data = np.linalg.norm(data, axis=2)

    ref = import_pestoseis_data(
        f"reference_solutions/pestoseis_elastic_hom_{suffix}_ref.h5"
    )["waveforms"]

    # Ensure sampling between the reference data and the simulated data are the same
    ref_resampled = np.zeros_like(data)
    for i, trace in enumerate(ref):
        interp_fun = interpolate.interp1d(np.linspace(0, 1, ref.shape[1]), trace)
        ref_resampled[i, :] = interp_fun(np.linspace(0, 1, data.shape[1]))

    diff = np.mean(np.linalg.norm(data - ref_resampled, axis=1))
    print(
        f"L2 Norm between simulated data and the PestoSeis reference solution is {diff}\n"
    )

    return diff < threshold


def _test_elastic_gradient(threshold: float = 5e-9) -> bool:
    """
    Test function for checking if the elastic wave solver produces the expected
    results.

    Solves a simple forward problem of a wave propagating through a a medium
    with linearly varying material properties with depth.  The result is
    compared to a reference solution within the `reference_solutions` directory.
    If the difference (measured in an L2 sense) is greater than the user-defined
    threshold, then the test will fail.

    Args:
        threshold (float, optional): Acceptable error threshold above which the
        test will fail.

    Returns:
        (bool): Whether the test passed or failed.
    """
    vp_model = np.zeros(domain_size, dtype=float)
    vs_model = np.zeros(domain_size, dtype=float)
    rho_model = np.zeros(domain_size, dtype=float)
    for i in range(nx):
        vp_model[i, :] = np.linspace(vp_ref_range[0], vp_ref_range[1], nz)
        vs_model[i, :] = np.linspace(vs_ref_range[0], vs_ref_range[1], nz)
        rho_model[i, :] = np.linspace(rho_ref_range[0], rho_ref_range[1], nz)

    lambda_model = rho_model * (vp_model ** 2 - 2.0 * vs_model ** 2)
    mu_model = rho_model * vs_model ** 2

    rockprops = {}
    rockprops["lambda"] = lambda_model
    rockprops["mu"] = mu_model
    rockprops["rho"] = rho_model

    inpar["freesurface"] = False
    suffix = "PML"

    inpar["sourcetype"] = "MomTensor"
    inpar["momtensor"] = dict(xx=1.0, xz=0.0, zz=1.0)
    inpar["seismogrkind"] = "velocity"

    stf = pestoseis.seismicwaves2d.sourcetimefuncs.ricker_1st_derivative_source(
        t, t0, f0
    )

    data, _ = pestoseis.seismicwaves2d.elasticwaveprop2D.solveelastic2D(
        inpar,
        rockprops,
        sources,
        stf,
        f0,
        receivers,
        outfileh5=f"pestoseis_elastic_grad_{suffix}.h5",
    )

    data = np.linalg.norm(data, axis=2)

    ref = import_pestoseis_data(
        f"reference_solutions/pestoseis_elastic_grad_{suffix}_ref.h5"
    )["waveforms"]

    # Ensure sampling between the reference data and the simulated data are the same
    ref_resampled = np.zeros_like(data)
    for i, trace in enumerate(ref):
        interp_fun = interpolate.interp1d(np.linspace(0, 1, ref.shape[1]), trace)
        ref_resampled[i, :] = interp_fun(np.linspace(0, 1, data.shape[1]))

    diff = np.mean(np.linalg.norm(data - ref_resampled, axis=1))
    print(
        f"L2 Norm between simulated data and the PestoSeis reference solution is {diff}\n"
    )

    return diff < threshold


class TestAcousticHomogeneous(unittest.TestCase):
    def test_PML(self):
        self.assertTrue(
            _test_acoustic_homogeneous(free_surface=False),
            "Test data does not match the saved reference data.",
        )

    def test_free_surface(self):
        self.assertTrue(
            _test_acoustic_homogeneous(free_surface=True),
            "Test data does not match the saved reference data.",
        )

    def test_PML_gradient(self):
        self.assertTrue(
            _test_acoustic_gradient(),
            "Test data does not match the saved reference data.",
        )


class TestElasticHomogeneous(unittest.TestCase):
    def test_PML(self):
        self.assertTrue(
            _test_elastic_homogeneous(free_surface=False),
            "Test data does not match the saved reference data.",
        )

    def test_free_surface(self):
        self.assertTrue(
            _test_elastic_homogeneous(free_surface=True),
            "Test data does not match the saved reference data.",
        )

    def test_PML_gradient(self):
        self.assertTrue(
            _test_elastic_gradient(),
            "Test data does not match the saved reference data.",
        )


# -

if __name__ == "__main__":
    unittest.main()

    # Functions for generating the test data; must be manually enabled to
    # re-compute the reference data
    if False:
        # Acoustic reference data
        _acoustic_ref_data(free_surface=False)
        _acoustic_ref_data(free_surface=True)
        _acoustic_ref_data(free_surface=False, homogeneous=False)

        # Elastic reference data
        _elastic_ref_data(free_surface=False)
        _elastic_ref_data(free_surface=True)
        _elastic_ref_data(free_surface=False, homogeneous=False)
