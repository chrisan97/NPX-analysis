# Written on 5/29/2024
# Chris (Seong Yeol) An
# san11@jh.edu

1. First run CatGT_runFile.
	- You will have to adjust the parameters accordingly to read the files that you want.
	- CatGT will do several things.
	1. Apply CAR over the recording file. 
	2. Extract recording edge files. This will essentially extract the important edges.

2. Sort the spikes using KS2.5 MATLAB version to the kilosort_output folder.
	- KS2.5 will extract the spikes, generate all the necessary .npy files for sorting.

3. Manually curate the spikes in phy, using the outputs from KS2.5/

4. After you sort the spikes in phy, take the spikes.npy (which has the spikes from KS2.5).
	- You now need to run TPrime_runFile_Jupyter to shift the time windows.

