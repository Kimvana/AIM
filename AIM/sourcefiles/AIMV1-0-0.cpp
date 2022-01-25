// AIM_V9.cpp : Defines the exported functions for the DLL application.
//

#include <math.h>
#include <stdio.h>

/*extern "C" {
	__declspec(dllexport) void CalcFieldGrad(float *refpos, int *inrange_ix, float *WS_positions, float *WS_charges, float *halfbox, float *boxdims, int Natoms, float *out);
	__declspec(dllexport) void CalcField(float *refpos, int *inrange_ix, float *WS_positions, float *WS_charges, float *halfbox, float *boxdims, int Natoms, float *out);
	__declspec(dllexport) void group_difference(int *atlist1, int *atlist2, int *lengths);
	__declspec(dllexport) void res_finder(int *COM_resnums, int COM_resnums_len, int *residues_start, int *residues_fin, float *WS_positions, float *WS_masses, float *WS_res_COM, float *halfbox, float *boxdims);
    __declspec(dllexport) void res_finder_dumb(int *COM_resnums, int COM_resnums_len, int *residues_start, int *residues_fin, float *WS_positions, float *WS_masses, float *WS_res_COM, float *halfbox, float *boxdims);
	__declspec(dllexport) int AGsorter(int refind, int *selection_residues_resnums, int sel_res_resnums_len, float *WS_res_COM, int *inrange_res, float desrange, int WS_nres, float *halfbox, float *boxdims);
	__declspec(dllexport) int ix_builder(int *inrange_reslist, int inrange_reslist_len, int *residues_start, int *residues_fin, int *inrange_ix);
	__declspec(dllexport) void PrepTCC(int *AllAmGroups, int *PrePro, float *TCC_v, int *relevant_groups, int relgrouplen, float *vPro, float *vGen, float alphaPro, float alphaGen, float *WS_positions, float *halfbox, float *boxdims);
	__declspec(dllexport) float CalcTCC(int AmGroupi, int AmGroupj, int *AllAmGroups, int *PrePro, float *TCC_v, float *qGen, float *qPro, float *dqGen, float *dqPro, float TCC_4PiEps, float *WS_positions, float *halfbox, float *boxdims);
	__declspec(dllexport) void PrepTDCKrimm(int *AllAmGroups, float *TDCKr_r, float *TDCKr_m, int *relevant_groups, int relgrouplen, float *WS_positions, float *halfbox, float *boxdims);
	__declspec(dllexport) float CalcTDCKrimm(int AmGroupi, int AmGroupj, float *TDCKr_r, float *TDCKr_m, float *halfbox, float *boxdims);
	__declspec(dllexport) void PrepTDCTasumi(int *AllAmGroups, float *TDCTa_r, float *TDCTa_m, int *relevant_groups, int relgrouplen, float *WS_positions, float *halfbox, float *boxdims);
	__declspec(dllexport) float CalcTDCTasumi(int AmGroupi, int AmGroupj, float *TDCTa_r, float *TDCTa_m, float *halfbox, float *boxdims);
    __declspec(dllexport) float CalcGenCoup(int OscGroupi, int OscGroupj, float *TDCGen_r, float *TDCGen_m, float *halfbox, float *boxdims);
}*/


extern "C" {
	// (for now) internal-only math functions. (if external, add to list above, at least for windows)
	void PBC_diff(float *vect1, float *vect2, float *halfbox, float *boxdims, float *vectout) {
		for (int i = 0; i < 3; i++) {
			vectout[i] = vect1[i] - vect2[i];
			//vectout[i] = (vectout[i] + halfbox[i]) % boxdims[i] - halfbox[i];
			if (vectout[i] > halfbox[i]) {
				vectout[i] -= boxdims[i];
			}
			else if (vectout[i] < -1 * halfbox[i]) {
				vectout[i] += boxdims[i];
			}
		}
	}

	float dotprod(float *vect1, float *vect2) {
		float out = vect1[0] * vect2[0] + vect1[1] * vect2[1] + vect1[2] * vect2[2];
		return out;
	}

	float veclen(float *vect) {
		return sqrt(dotprod(vect, vect));
	}

	void vecnorm(float *vect) {
		float iveclen = 1/veclen(vect);
		for (int i = 0; i < 3; i++) {
			vect[i] *= iveclen;
		}
	}

	void project(float *vect1, float *vect2) {
		float inprod = dotprod(vect1, vect2) / dotprod(vect1, vect1);
		for (int i = 0; i < 3; i++) {
			vect2[i] -= inprod * vect1[i];
		}
	}

	void crossprod(float *vect1, float *vect2, float *vectout) {
		vectout[0] = vect1[1] * vect2[2] - vect1[2] * vect2[1];
		vectout[1] = vect1[2] * vect2[0] - vect1[0] * vect2[2];
		vectout[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];
	}


	void CalcFieldGrad(float *refpos, int *inrange_ix, float *WS_positions, float *WS_charges, float *halfbox, float *boxdims, int Natoms, float *out) {
		for (int i = 0; i < 10; i++) {
			out[i] = 0.0;
		}
		int dataloc;
		float charge, dist2, idist2, idist, prefac, prefac2, diffX, diffY, diffZ;
		for (int atom = 0; atom < Natoms; atom++) {
			float diff[3];
			dataloc = inrange_ix[atom];
			if (dataloc == 999999999) break;
			PBC_diff(refpos, &WS_positions[dataloc * 3], halfbox, boxdims, diff);
			charge = WS_charges[dataloc];
			dist2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
			idist2 = 1 / dist2;
			idist = sqrt(idist2);
			prefac = charge * idist2 * idist;
			prefac2 = 3.0 * prefac * idist2;

			diffX = diff[0];
			diffY = diff[1];
			diffZ = diff[2];

			out[0] += charge * idist;

			out[1] += diffX * prefac;
			out[2] += diffY * prefac;
			out[3] += diffZ * prefac;

			out[4] += prefac - (diffX * diffX * prefac2);
			out[5] += prefac - (diffY * diffY * prefac2);
			out[6] += prefac - (diffZ * diffZ * prefac2);
			out[7] -= diffX * diffY * prefac2;
			out[8] -= diffX * diffZ * prefac2;
			out[9] -= diffY * diffZ * prefac2;
		}
	}

	void CalcField(float *refpos, int *inrange_ix, float *WS_positions, float *WS_charges, float *halfbox, float *boxdims, int Natoms, float *out) {
		for (int i = 0; i < 4; i++) {
			out[i] = 0.0;
		}
		for (int atom = 0; atom < Natoms; atom++) {
			float diff[3];
			int dataloc = inrange_ix[atom];
			if (dataloc == 999999999) break;
			PBC_diff(refpos, &WS_positions[dataloc * 3], halfbox, boxdims, diff);

			float charge = WS_charges[dataloc];
			float dist2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
			float idist2 = 1 / dist2;
			float idist = sqrt(idist2);
			float prefac = charge * idist2 * idist;

			float diffX = diff[0];
			float diffY = diff[1];
			float diffZ = diff[2];

			out[0] += charge * idist;

			out[1] += diffX * prefac;
			out[2] += diffY * prefac;
			out[3] += diffZ * prefac;
		}
	}

	void group_difference(int *atlist1, int *atlist2, int *lengths) {
		int Nat1 = lengths[0];
		int Nat2 = lengths[1];

		int at1c = 0;
		int at2c = 0;
		int at3c = 0;

		while (true) {
			if (atlist1[at1c] > atlist2[at2c]) {
				at2c++;
			}
			else if (atlist1[at1c] < atlist2[at2c]) {
				atlist1[at3c] = atlist1[at1c];
				at3c++;
				at1c++;
			}
			else if (atlist1[at1c] == atlist2[at2c]) {
				at1c++;
				at2c++;
			}

			if (at1c == Nat1) {
				break;
			}
			else if (at2c == Nat2) {
				for (int i = at1c; i < Nat1; i++) {
					atlist1[at3c] = atlist1[i];
					at3c++;
				}

				break;
			}
		}

		for (int i = at3c; i < Nat1; i++) {
			atlist1[i] = 999999999;
		}

		lengths[2] = at3c;
	}

	void res_finder(int *COM_resnums, int COM_resnums_len, int *residues_start, int *residues_fin, float *WS_positions, float *WS_masses, float *WS_res_COM, float *halfbox, float *boxdims) {
		for (int index = 0; index < COM_resnums_len; index++) {
			int resnum = COM_resnums[index];
			int start = residues_start[resnum];
			int fin = residues_fin[resnum];

			float refpos[3];
			float cumpos[3];
			for (int i = 0; i < 3; i++) {
				refpos[i] = WS_positions[start * 3 + i];
				cumpos[i] = 0.0;
			}
			float cummass = WS_masses[start];

			for (int atom = (start + 1); atom < (fin + 1); atom++) {
				float diff[3];
				float atmass = WS_masses[atom];
				cummass += atmass;
				for (int i = 0; i < 3; i++) {
					diff[i] = WS_positions[atom * 3 + i] - refpos[i];
					if (diff[i] > halfbox[i]) {
						diff[i] -= boxdims[i];
					}
					else if (diff[i] < -1 * halfbox[i]) {
						diff[i] += boxdims[i];
					}
					cumpos[i] += (diff[i] * atmass);
				}
			}

			float icummass = 1 / cummass;

			/*for (int i = 0; i < 3; i++) {
				cumpos[i] = (cumpos[i] * icummass) + refpos[i];
				if (cumpos[i] > halfbox[i]) {
					cumpos[i] -= boxdims[i];
				}
				else if (cumpos[i] < -1 * halfbox[i]) {
					cumpos[i] += boxdims[i];
				}
				WS_res_COM[resnum * 3 + i] = cumpos[i];
			}*/

			for (int i = 0; i < 3; i++) {
				cumpos[i] = (cumpos[i] * icummass) + refpos[i];
				if (cumpos[i] > boxdims[i]) {
					cumpos[i] -= boxdims[i];
				}
				else if (cumpos[i] < 0) {
					cumpos[i] += boxdims[i];
				}
				WS_res_COM[resnum * 3 + i] = cumpos[i];
			}

		}
	}

	void res_finder_dumb(int *COM_resnums, int COM_resnums_len, int *residues_start, int *residues_fin, float *WS_positions, float *WS_masses, float *WS_res_COM, float *halfbox, float *boxdims) {
		for (int index = 0; index < COM_resnums_len; index++) {
			int resnum = COM_resnums[index];
			int start = residues_start[resnum];
			int fin = residues_fin[resnum];

			float cumpos[3];
			for (int i = 0; i < 3; i++) {
				cumpos[i] = 0.0;
			}
			float cummass = 0.0;

			for (int atom = (start); atom < (fin + 1); atom++) {
				float diff[3];
				float atmass = WS_masses[atom];
				cummass += atmass;
				for (int i = 0; i < 3; i++) {
					diff[i] = WS_positions[atom * 3 + i];
					cumpos[i] += (diff[i] * atmass);
				}
			}

			float icummass = 1 / cummass;

			for (int i = 0; i < 3; i++) {
				cumpos[i] = (cumpos[i] * icummass);
				WS_res_COM[resnum * 3 + i] = cumpos[i];
			}

		}
	}

	int AGsorter(int refind, int *selection_residues_resnums, int sel_res_resnums_len, float *WS_res_COM, int *inrange_res, float desrange, int WS_nres, float *halfbox, float *boxdims) {
		float desrange2 = desrange * desrange;
		float refCOM[3];
		int counter = 0;
		for (int i = 0; i < 3; i++) {
			refCOM[i] = WS_res_COM[refind * 3 + i];
		}
		for (int index = 0; index < sel_res_resnums_len; ++index) {
			int resnum = selection_residues_resnums[index];
			float diff[3];
			PBC_diff(refCOM, &WS_res_COM[resnum * 3], halfbox, boxdims, diff);

			float dist2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

			if (dist2 < desrange2) {
				inrange_res[counter] = resnum;
				counter += 1;
			}
		}
		for (int i = counter; i < WS_nres; i++) {
			inrange_res[i] = 999999999;
		}
		return counter;
	}

	int ix_builder(int *inrange_reslist, int inrange_reslist_len, int *residues_start, int *residues_fin, int *inrange_ix) {
		int counter = 0;
		for (int index = 0; index < inrange_reslist_len; index++) {
			int resnum = inrange_reslist[index];
			int start = residues_start[resnum];
			int fin = residues_fin[resnum];
			for (int atnum = start; atnum < (fin + 1); atnum++) {
				inrange_ix[counter] = atnum;
				counter++;
			}
		}
		return counter;
	}

	void PrepTCC(int *AllAmGroups, int *PrePro, float *TCC_v, int *relevant_groups, int relgrouplen, float *vPro, float *vGen, float alphaPro, float alphaGen, float *WS_positions, float *halfbox, float *boxdims) {
		int C, O, N, AmGroup;
		float CO[3], CN[3], z[3];
        for (int count = 0; count < relgrouplen; count++) {
            AmGroup = relevant_groups[count];
			C = AllAmGroups[AmGroup * 18 + 6];
			O = AllAmGroups[AmGroup * 18 + 7];
			N = AllAmGroups[AmGroup * 18 + 9];
			PBC_diff(&WS_positions[O * 3], &WS_positions[C * 3], halfbox, boxdims, CO);
			PBC_diff(&WS_positions[N * 3], &WS_positions[C * 3], halfbox, boxdims, CN);
			vecnorm(CO);
			project(CO, CN);
			vecnorm(CN);
			crossprod(CO, CN, z);

			if (PrePro[AmGroup] == 0) {
				for (int a = 0; a < 6; a++) {
					for (int c = 0; c < 3; c++) {
						TCC_v[AmGroup * 18 + a * 3 + c] = (CO[c] * vGen[a * 3] + CN[c] * vGen[a * 3 + 1] + z[c] * vGen[a * 3 + 2])*alphaGen;
					}
				}
			}
			else {
				for (int a = 0; a < 6; a++) {
					for (int c = 0; c < 3; c++) {
						TCC_v[AmGroup * 18 + a * 3 + c] = (CO[c] * vPro[a * 3] + CN[c] * vPro[a * 3 + 1] + z[c] * vPro[a * 3 + 2])*alphaPro;
					}
				}
			}
		}

	}

	float CalcTCC(int AmGroupi, int AmGroupj, int *AllAmGroups, int *PrePro, float *TCC_v, float *qGen, float *qPro, float *dqGen, float *dqPro, float TCC_4PiEps, float *WS_positions, float *halfbox, float *boxdims) {
		float J = 0;
		int atom_a, atom_b;
		float va[3], vb[3], d[3];
		float r2, ir2, ir, ir3, ir5;
		float qi, qj[6], dqi, dqj[6];
		int tyi, tyj;

		tyi = PrePro[AmGroupi];
		tyj = PrePro[AmGroupj];

		for (int b = 0; b < 6; b++) {
			qj[b] = (tyj * qPro[b] + (1-tyj) * qGen[b]);
			dqj[b] = (tyj * dqPro[b] + (1-tyj) * dqGen[b]);
		}

		for (int a = 0; a < 6; a++) {
			qi = (tyi * qPro[a] + (1-tyi) * qGen[a]);
			dqi = (tyi * dqPro[a] + (1-tyi) * dqGen[a]);
			atom_a = AllAmGroups[AmGroupi * 18 + 6 + a];
			for (int c = 0; c < 3; c++) {
				va[c] = TCC_v[AmGroupi * 18 + a * 3 + c];
			}
			for (int b = 0; b < 6; b++) {
				atom_b = AllAmGroups[AmGroupj * 18 + 6 + b];
				PBC_diff(&WS_positions[atom_b * 3], &WS_positions[atom_a * 3], halfbox, boxdims, d);
				r2 = dotprod(d, d);
				if (r2 < 0.01) r2 = 1.00;
				ir2 = 1 / r2;
				ir = sqrt(ir2);
				ir3 = ir * ir2;
				ir5 = ir3 * ir2;

				for (int c = 0; c < 3; c++) {
					vb[c] = TCC_v[AmGroupj * 18 + b * 3 + c];
				}

				J -= 3 * ir5 * qi * qj[b] * dotprod(vb, d) * dotprod(va, d);
				J -= ir3 * (dqi * qj[b] * dotprod(vb, d) - qi * dqj[b] * dotprod(va, d) - dotprod(va, vb) * qi * qj[b]);
				J += ir * dqi * dqj[b];

			}
		}
		J *= TCC_4PiEps;
		return J;
	}

	void PrepTDCKrimm(int *AllAmGroups, float *TDCKr_r, float *TDCKr_m, int *relevant_groups, int relgrouplen, float *WS_positions, float *halfbox, float *boxdims) {
		int C, O, N, AmGroup;
		float CO[3], CN[3], m[3], dispibond, bond;
		float displace = 0.868;
		float tanangle = -0.36397023426;   // = -tan(20.0 / 180.0*3.14159265359);
        
        for (int count = 0; count < relgrouplen; count++) {
            AmGroup = relevant_groups[count];
			C = AllAmGroups[AmGroup * 18 + 6];
			O = AllAmGroups[AmGroup * 18 + 7];
			N = AllAmGroups[AmGroup * 18 + 9];
			PBC_diff(&WS_positions[O * 3], &WS_positions[C * 3], halfbox, boxdims, CO);
			dispibond = displace/veclen(CO);
			for (int i = 0; i < 3; i++) {
				TDCKr_r[AmGroup * 3 + i] = WS_positions[C * 3 + i] + dispibond * CO[i];
			}
			PBC_diff(&WS_positions[N * 3], &WS_positions[C * 3], halfbox, boxdims, CN);
			vecnorm(CO);
			project(CO, CN);
			vecnorm(CN);
			bond = tanangle;
			for (int i = 0; i < 3; i++) {
				m[i] = CO[i] + CN[i] * bond;
			}
			vecnorm(m);
			for (int i = 0; i < 3; i++) {
				TDCKr_m[AmGroup * 3 + i] = m[i];
			}
		}
	}

	float CalcTDCKrimm(int AmGroupi, int AmGroupj, float *TDCKr_r, float *TDCKr_m, float *halfbox, float *boxdims) {
		float J = 0;
		//float fourPiEps = 5034.13; //219476*0.529177249*0.208194*0.208194
		float fourPiEps = 580;
		float d[3];
		PBC_diff(&TDCKr_r[AmGroupi * 3], &TDCKr_r[AmGroupj * 3], halfbox, boxdims, d);
		float ir = 1/veclen(d);
		float ir2 = ir * ir;
		float ir3 = ir * ir2;
		float ir5 = ir3 * ir2;
		J = fourPiEps * (dotprod(&TDCKr_m[AmGroupi * 3], &TDCKr_m[AmGroupj * 3]) * ir3 - 3.0 * dotprod(&TDCKr_m[AmGroupi * 3], d) * dotprod(&TDCKr_m[AmGroupj * 3], d) * ir5);
		return J;
	}

	void PrepTDCTasumi(int *AllAmGroups, float *TDCTa_r, float *TDCTa_m, int *relevant_groups, int relgrouplen, float *WS_positions, float *halfbox, float *boxdims) {
		int C, O, N, AmGroup;
		float CO[3], CN[3], dri[3], COvecDri, prefac, m[3];
		float itheta = 5.6712818196; // = 1/tan(0.17632698 = 10 deg expressed in rad)

		// float magnitude = 2.73 * sqrt(51.43/5034);
        float magnitude = 0.276;

        for (int count = 0; count < relgrouplen; count++) {
            AmGroup = relevant_groups[count];
			C = AllAmGroups[AmGroup * 18 + 6];
			O = AllAmGroups[AmGroup * 18 + 7];
			N = AllAmGroups[AmGroup * 18 + 9];
			PBC_diff(&WS_positions[O * 3], &WS_positions[C * 3], halfbox, boxdims, CO);
			PBC_diff(&WS_positions[N * 3], &WS_positions[C * 3], halfbox, boxdims, CN);
			vecnorm(CO);
			vecnorm(CN);

			for (int i = 0; i < 3; i++) {
				dri[i] = 0.665*CO[i] + 0.258 * CN[i];
				TDCTa_r[AmGroup * 3 + i] = WS_positions[C * 3 + i] + dri[i];
			}

			COvecDri = dotprod(CO, dri);
			prefac = COvecDri + sqrt(dotprod(dri, dri) - COvecDri * COvecDri)*itheta;
			for (int i = 0; i < 3; i++) {
				m[i] = dri[i] - prefac * CO[i];
			}

			vecnorm(m);
			for (int i = 0; i < 3; i++) {
				TDCTa_m[AmGroup * 3 + i] = m[i] * magnitude;
			}
		}
	}

	float CalcTDCTasumi(int AmGroupi, int AmGroupj, float *TDCTa_r, float *TDCTa_m, float *halfbox, float *boxdims) {
		float J = 0;
		// float AoverEps = 51.43; // =0.1*848619/1650
        float AoverEps = 5034;  // see comment in couplingfunctions.py - CalcGenCoup
		float d[3];
		PBC_diff(&TDCTa_r[AmGroupi * 3], &TDCTa_r[AmGroupj * 3], halfbox, boxdims, d);
		float ir = 1 / veclen(d);
		float ir2 = ir * ir;
		float ir3 = ir * ir2;
		float ir5 = ir3 * ir2;
		J = AoverEps * (dotprod(&TDCTa_m[AmGroupi * 3], &TDCTa_m[AmGroupj * 3]) * ir3 - 3.0 * dotprod(&TDCTa_m[AmGroupi * 3], d) * dotprod(&TDCTa_m[AmGroupj * 3], d) * ir5);
		return J;
	}

    float CalcGenCoup(int OscGroupi, int OscGroupj, float *TDCGen_r, float *TDCGen_m, float *halfbox, float *boxdims){
        float J = 0;
        float fourPiEps = 5034;  // see CoupFuncs.py - CalcGenCoup_nb for explanation!
        float d[3];
        PBC_diff(&TDCGen_r[OscGroupi * 3], &TDCGen_r[OscGroupj * 3], halfbox, boxdims, d);
        float ir = 1/veclen(d);
        float ir2 = ir * ir;
        float ir3 = ir * ir2;
        float ir5 = ir3 * ir2;
        J = fourPiEps * (dotprod(&TDCGen_m[OscGroupi * 3], &TDCGen_m[OscGroupj * 3]) * ir3 - 3.0 * dotprod(&TDCGen_m[OscGroupi * 3], d) * dotprod(&TDCGen_m[OscGroupj * 3], d) * ir5);
        return J;
    }

}



