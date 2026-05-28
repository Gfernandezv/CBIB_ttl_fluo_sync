import pyabf
import matplotlib.pyplot as plt
import pyabf.tools.memtest
import numpy as np
import json

def read_abf(file):
    return pyabf.ABF(file)

def compute_memtest(abf):
    memtest = pyabf.tools.memtest.Memtest(abf)
    return {
        "sweepTimesMin": abf.sweepTimesMin.tolist(),
        "Ih": memtest.Ih.values.tolist(),
        "Rm": memtest.Rm.values.tolist(),
        "Ra": memtest.Ra.values.tolist(),
        "CmStep": memtest.CmStep.values.tolist()
    }

def extract_patch_data(abf, pt1=259, pt2=259.5, pt3=255, pt4=257):
    currents, voltages = [], []
    for sweep in abf.sweepList:
        abf.setSweep(sweep, channel=2)
        currents.append(np.mean(abf.sweepY[int(pt1 * abf.dataPointsPerMs): int(pt2 * abf.dataPointsPerMs)]) * -1)
        abf.setSweep(sweep, channel=1)
        voltages.append(np.mean(abf.sweepY[int(pt3 * abf.dataPointsPerMs): int(pt4 * abf.dataPointsPerMs)]))
    return np.array(currents), np.array(voltages)

def plot_traces(abf, pt1, pt2, pt3, pt4):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    
    for sweep in abf.sweepList:
        abf.setSweep(sweep, channel=2)
        ax1.plot(abf.sweepX, abf.sweepY, lw=0.5)
    ax1.set_ylabel(abf.sweepLabelY)
    
    for sweep in abf.sweepList:
        abf.setSweep(sweep, channel=1)
        ax2.plot(abf.sweepX, abf.sweepY, lw=0.5)
    ax2.set_xlabel(abf.sweepLabelX)
    ax2.set_ylabel(abf.sweepLabelC)
    
    for ax in (ax1, ax2):
        for pt in [pt1, pt2, pt3, pt4]:
            ax.axvline(pt / 10000, alpha=0.5, color='k', ls='--', lw=0.5)

    plt.tight_layout()
    plt.show()

def analize_patch(file, graph=True):
    abf = read_abf(file)
    memtest_data = compute_memtest(abf)
    
    currents, voltages = extract_patch_data(abf)
    corrected_currents = currents / np.array(memtest_data['CmStep'])
    
    results = {
        "patch_data": {
            "currents": currents.tolist(),
            "voltages": voltages.tolist(),
            "corrected_currents": corrected_currents.tolist()
        },
        "memtest_data": memtest_data
    }
    
    if graph:
        plot_traces(abf, 259, 259.5, 255, 257)
    
    # Guardar resultados opcionalmente
    with open(file.replace(".abf", "_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    
    return results

from pynwb import NWBFile, NWBHDF5IO
from pynwb.file import Subject
from pynwb.icephys import VoltageClampSeries, VoltageClampStimulusSeries
import datetime
import numpy as np
import os

def save_to_nwb_voltageclamp_full(file, output_path):
    abf = pyabf.ABF(file)

    nwbfile = NWBFile(
        session_description='Whole-cell voltage clamp in HEK cells',
        identifier=os.path.basename(file),
        session_start_time=datetime.datetime.now().astimezone(),
    )

    device = nwbfile.create_device(name='PatchClampRig')
    electrode = nwbfile.create_icephys_electrode(
        name="electrode_1",
        description="Whole-cell voltage clamp electrode",
        device=device
    )

    for sweep in abf.sweepList:
        # ⚡ Corriente cruda (Channel 0)
        abf.setSweep(sweep, channel=0)
        current_raw = abf.sweepY.copy()
        current_series = VoltageClampSeries(
            name=f"recorded_current_raw_sweep_{sweep}",
            data=current_raw,
            unit='pA',
            electrode=electrode,
            conversion=1.0,
            resolution=np.nan,
            starting_time=float(abf.sweepX[0]),
            rate=float(abf.dataRate)
        )
        nwbfile.add_acquisition(current_series)

        # 🔌 Voltaje de comando (Channel 1)
        abf.setSweep(sweep, channel=1)
        voltage_command = abf.sweepY.copy()
        voltage_command_series = VoltageClampStimulusSeries(
            name=f"commanded_voltage_sweep_{sweep}",
            data=voltage_command,
            unit='mV',
            electrode=electrode,
            conversion=1.0,
            resolution=np.nan,
            starting_time=float(abf.sweepX[0]),
            rate=float(abf.dataRate)
        )
        nwbfile.add_stimulus(voltage_command_series)

        # 📊 Canal 2 (corriente leak-subtracted)
        abf.setSweep(sweep, channel=2)
        leak_corrected = abf.sweepY.copy()
        aux_series = VoltageClampSeries(
            name=f"recorded_current_leak_subtracted_sweep_{sweep}",
            data=leak_corrected,
            unit='pA',
            electrode=electrode,
            conversion=1.0,
            resolution=np.nan,
            starting_time=float(abf.sweepX[0]),
            rate=float(abf.dataRate),
            comments="Leak-subtracted current"
        )
        nwbfile.add_acquisition(aux_series)

    output_file = os.path.join(output_path, os.path.basename(file).replace(".abf", ".nwb"))
    with NWBHDF5IO(output_file, 'w') as io:
        io.write(nwbfile)

    print(f"NWB file saved with raw + voltage command + leak-subtracted sweeps: {output_file}")
