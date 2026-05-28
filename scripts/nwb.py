import os
import glob
import json
from datetime import datetime
from dateutil.tz import tzlocal
import numpy as np
import pyabf
import pyabf.tools.memtest
import pandas as pd

from pynwb import NWBHDF5IO, NWBFile
from pynwb.file import Subject
from pynwb.icephys import VoltageClampStimulusSeries, VoltageClampSeries, IntracellularRecordingsTable
from pynwb.misc import DynamicTable
from hdmf.backends.hdf5.h5_utils import H5DataIO


DEFAULT_METADATA = {
    "experimenter": ["Germán Fernández"],
    "experiment_description": "Whole-cell patch-clamp recording in HEK293 cells under TRPM3 protocol.",
    "institution": "Universidad Nacional Andrés Bello",
    "keywords": ["patch-clamp", "voltage clamp", "HEK293", "TRPM3"],
    "subject_id": "HEK293_01",
    "subject_description": "HEK293 cell expressing TRPM3 channels",
    "device_description": "Multiclamp 700B amplifier with HEKA EPC, 10 kHz sampling",
    "electrode_description": "Whole-cell voltage clamp electrode",
    "cell_id": "HEK293_Cell_001"
}


def createCompressedDataset(array):
    """Request compression for the given array and return it wrapped."""
    return H5DataIO(data=array, compression=True, chunks=True, shuffle=True, fletcher32=True)


class ABF2Converter:
    def __init__(self, inputPath, outputFilePath, include_leak=True, include_memtest=True, metadata=None):
        self.inputPath = inputPath
        self.include_leak = include_leak
        self.include_memtest = include_memtest
        self.metadata = metadata or DEFAULT_METADATA

        if os.path.isfile(self.inputPath):
            abf = pyabf.ABF(self.inputPath)
            if abf.abfVersion["major"] != 2:
                raise ValueError("Only ABF v2.x supported.")
            self.fileNames = [os.path.basename(self.inputPath)]
            self.abfFiles = [abf]

        elif os.path.isdir(self.inputPath):
            abfFiles = sorted(glob.glob(os.path.join(self.inputPath, "*.abf")))
            if len(abfFiles) == 0:
                raise ValueError(f"No ABF files found in {self.inputPath}.")
            self.fileNames = [os.path.basename(f) for f in abfFiles]
            self.abfFiles = [pyabf.ABF(f) for f in abfFiles]
        else:
            raise ValueError("Invalid input path.")

        self.outputPath = outputFilePath

    def _createNWBFile(self):
        start_time = self.abfFiles[0].abfDateTime
        self.NWBFile = NWBFile(
            session_description=self.metadata["experiment_description"],
            session_start_time=start_time,
            identifier="ABF2Session",
            file_create_date=datetime.now(tzlocal()),
            experimenter=self.metadata["experimenter"],
            experiment_description=self.metadata["experiment_description"],
            institution=self.metadata["institution"],
            keywords=self.metadata["keywords"]
        )
        self.NWBFile.subject = Subject(
            subject_id=self.metadata["subject_id"],
            description=self.metadata["subject_description"]
        )
        return self.NWBFile

    def _createDeviceAndElectrode(self):
        self.device = self.NWBFile.create_device(
            name="PatchClampRig",
            description=self.metadata["device_description"]
        )
        self.electrode = self.NWBFile.create_icephys_electrode(
            name='elec0',
            device=self.device,
            description=self.metadata["electrode_description"],
            cell_id=self.metadata["cell_id"]
        )

    def _unitConversion(self, unit):
        """Return conversion factor and SI unit."""
        if unit.lower() == 'pa':
            return 1e-12, 'amperes'
        elif unit.lower() == 'mv':
            return 1e-3, 'volts'
        elif unit.lower() == 'v':
            return 1.0, 'volts'
        elif unit.lower() == 'a':
            return 1.0, 'amperes'
        return 1.0, unit

    def _addAcquisitionWithTable(self):
        """Add sweeps independently and register in IntracellularRecordingsTable."""
        recordings_table = IntracellularRecordingsTable()

        for idx, abf in enumerate(self.abfFiles):
            for sweep in abf.sweepList:
                # Stimulus (Ch1)
                abf.setSweep(sweep, channel=1)
                conv_stim, unit_stim = self._unitConversion(abf.sweepUnitsY)
                stim_series = VoltageClampStimulusSeries(
                    name=f"stimulus_{idx}_{sweep}",
                    data=createCompressedDataset(abf.sweepY * conv_stim),
                    sweep_number=int(sweep),
                    electrode=self.electrode,
                    gain=1.0,
                    resolution=np.nan,
                    conversion=1.0,
                    starting_time=0.0,
                    rate=float(abf.dataRate),
                    unit=unit_stim
                )
                self.NWBFile.add_stimulus(stim_series)

                # Response (raw current, Ch0)
                abf.setSweep(sweep, channel=0)
                conv_resp, unit_resp = self._unitConversion(abf.sweepUnitsY)
                resp_series = VoltageClampSeries(
                    name=f"raw_current_{idx}_{sweep}",
                    data=createCompressedDataset(abf.sweepY * conv_resp),
                    sweep_number=int(sweep),
                    electrode=self.electrode,
                    gain=1.0,
                    resolution=np.nan,
                    conversion=1.0,
                    starting_time=0.0,
                    rate=float(abf.dataRate),
                    unit=unit_resp,
                    capacitance_fast=np.nan,
                    capacitance_slow=np.nan,
                    resistance_comp_bandwidth=np.nan,
                    resistance_comp_correction=np.nan,
                    resistance_comp_prediction=np.nan,
                    whole_cell_capacitance_comp=np.nan,
                    whole_cell_series_resistance_comp=np.nan
                )
                self.NWBFile.add_acquisition(resp_series)

                # Leak-subtracted current (optional)
                if self.include_leak and abf.channelCount > 2:
                    abf.setSweep(sweep, channel=2)
                    conv_leak, unit_leak = self._unitConversion(abf.sweepUnitsY)
                    leak_series = VoltageClampSeries(
                        name=f"leak_subtracted_{idx}_{sweep}",
                        data=createCompressedDataset(abf.sweepY * conv_leak),
                        sweep_number=int(sweep),
                        electrode=self.electrode,
                        gain=1.0,
                        resolution=np.nan,
                        conversion=1.0,
                        starting_time=0.0,
                        rate=float(abf.dataRate),
                        unit=unit_leak,
                        description="Leak-subtracted current"
                    )
                    self.NWBFile.add_acquisition(leak_series)

                # Register in intracellular_recordings table
                recordings_table.add_recording(
                    electrode=self.electrode,
                    stimulus=stim_series,
                    response=resp_series
                )

        self.NWBFile.intracellular_recordings = recordings_table

    def _addMemtest(self):
        """Add Memtest metrics as a DynamicTable in NWB analysis."""
        all_rows = []
        for idx, abf in enumerate(self.abfFiles):
            memtest = pyabf.tools.memtest.Memtest(abf)
            df = pd.DataFrame({
                "file": self.fileNames[idx],
                "sweep_time_min": abf.sweepTimesMin,
                "Ih": memtest.Ih.values,
                "Rm": memtest.Rm.values,
                "Ra": memtest.Ra.values,
                "CmStep": memtest.CmStep.values
            })
            all_rows.append(df)
        df_all = pd.concat(all_rows, ignore_index=True)

        table = DynamicTable(name='memtest_results', description='Memtest metrics per sweep')
        for col in df_all.columns:
            table.add_column(name=col, description=f"{col} values")
        for _, row in df_all.iterrows():
            table.add_row(**row.to_dict())

        analysis_module = self.NWBFile.create_processing_module(
            name="icephys",
            description="Intracellular electrophysiology analysis"
        )
        analysis_module.add(table)

    def convert(self):
        self._createNWBFile()
        self._createDeviceAndElectrode()
        self._addAcquisitionWithTable()
        if self.include_memtest:
            self._addMemtest()

        with NWBHDF5IO(self.outputPath, "w") as io:
            io.write(self.NWBFile, cache_spec=True)

        print(f"✅ Successfully converted to {self.outputPath}")
