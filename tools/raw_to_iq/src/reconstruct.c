#include <stdio.h>
#include <stdlib.h>
#include "packet_decode.h"

// p.78
float BRC0[4] = { 3.f,3.f,3.16f,3.53f };
float BRC1[4] = { 4.f,4.f,4.08f,4.37f };
float BRC2[6] = { 6.f,6.f,6.f,6.15f, 6.5f,6.88f };
float BRC3[7] = { 9.f,9.f,9.f,9.f,9.36f,9.50f, 10.1f };
float BRC4[9] = { 15.f,15.f,15.f,15.f,15.f,15.f, 15.22f, 15.50f, 16.05f };

// p.79
float NRL0[4] = { .3637f,1.0915f,1.8208f,2.6406f };
float NRL1[5] = { .3042f,.9127f,1.5216f,2.1313f,2.8426f };
float NRL2[7] = { .2305f,.6916f,1.1528f,1.6140f,2.0754f,2.5369f,3.1191f };
float NRL3[10] = { .1702f,.5107f,.8511f,1.1916f,1.5321f,1.8726f,2.2131f,2.5536f,2.8942f,3.3744f };
float NRL4[16] = { .1130f,.3389f,.5649f,.7908f,1.0167f,1.2428f,1.4687f,1.6947f,1.9206f,2.1466f,2.3725f,2.5985f,2.8244f,3.0504f,3.2764f,3.6623f };

// p.80
float SF[256] = { 0.f, 0.630f, 1.250f, 1.880f, 2.510f, 3.130f, 3.760f, 4.390f, 5.010f, 5.640f, 6.270f, 6.890f, 7.520f, 8.150f, 8.770f, 9.40f, 10.030f, 10.650f, 11.280f, 11.910f, 12.530f, 13.160f, 13.790f, 14.410f, 15.040f, 15.670f, 16.290f, 16.920f, 17.550f, 18.170f, 18.80f, 19.430f, 20.050f, 20.680f, 21.310f, 21.930f, 22.560f, 23.190f, 23.810f, 24.440f, 25.070f, 25.690f, 26.320f, 26.950f, 27.570f, 28.20f, 28.830f, 29.450f, 30.080f, 30.710f, 31.330f, 31.960f, 32.590f, 33.210f, 33.840f, 34.470f, 35.090f, 35.720f, 36.350f, 36.970f, 37.60f, 38.230f, 38.850f, 39.480f, 40.110f, 40.730f, 41.360f, 41.990f, 42.610f, 43.240f, 43.870f, 44.490f, 45.120f, 45.750f, 46.370f, 47.f, 47.630f, 48.250f, 48.880f, 49.510f, 50.130f, 50.760f, 51.390f, 52.010f, 52.640f, 53.270f, 53.890f, 54.520f, 55.150f, 55.770f, 56.40f, 57.030f, 57.650f, 58.280f, 58.910f, 59.530f, 60.160f, 60.790f, 61.410f, 62.040f, 62.980f, 64.240f, 65.490f, 66.740f, 68.f, 69.250f, 70.50f, 71.760f, 73.010f, 74.260f, 75.520f, 76.770f, 78.020f, 79.280f, 80.530f, 81.780f, 83.040f, 84.290f, 85.540f, 86.80f, 88.050f, 89.30f, 90.560f, 91.810f, 93.060f, 94.320f, 95.570f, 96.820f, 98.080f, 99.330f, 100.580f, 101.840f, 103.090f, 104.340f, 105.60f, 106.850f, 108.10f, 109.350f, 110.610f, 111.860f, 113.110f, 114.370f, 115.620f, 116.870f, 118.130f, 119.380f, 120.630f, 121.890f, 123.140f, 124.390f, 125.650f, 126.90f, 128.150f, 129.410f, 130.660f, 131.910f, 133.170f, 134.420f, 135.670f, 136.930f, 138.180f, 139.430f, 140.690f, 141.940f, 143.190f, 144.450f, 145.70f, 146.950f, 148.210f, 149.460f, 150.710f, 151.970f, 153.220f, 154.470f, 155.730f, 156.980f, 158.230f, 159.490f, 160.740f, 161.990f, 163.250f, 164.50f, 165.750f, 167.010f, 168.260f, 169.510f, 170.770f, 172.020f, 173.270f, 174.530f, 175.780f, 177.030f, 178.290f, 179.540f, 180.790f, 182.050f, 183.30f, 184.550f, 185.810f, 187.060f, 188.310f, 189.570f, 190.820f, 192.070f, 193.330f, 194.580f, 195.830f, 197.090f, 198.340f, 199.590f, 200.850f, 202.10f, 203.350f, 204.610f, 205.860f, 207.110f, 208.370f, 209.620f, 210.870f, 212.130f, 213.380f, 214.630f, 215.890f, 217.140f, 218.390f, 219.650f, 220.90f, 222.150f, 223.410f, 224.660f, 225.910f, 227.170f, 228.420f, 229.670f, 230.930f, 232.180f, 233.430f, 234.690f, 235.940f, 237.190f, 238.450f, 239.70f, 240.950f, 242.210f, 243.460f, 244.710f, 245.970f, 247.220f, 248.470f, 249.730f, 250.980f, 252.230f, 253.490f, 254.740f, 255.990f, 255.990f };

void reconstruction(unsigned char* BRCn, unsigned char* THIDXn, struct sh_code* hcode, int NQ, float* result)
{
	int hcode_index = 0, h;
	int BRCindex = 0;
	int inc = 128;
	do
	{
		if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
		for (h = 0; h < inc; h++) // p.68: 128 HCodes
		{
			switch (BRCn[BRCindex])
			{
			case 0:
				if (THIDXn[BRCindex] <= 3)
				{
					if (hcode[hcode_index].mcode < 3)
						result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
					else
						result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC0[THIDXn[BRCindex]];
				}
				else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL0[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
				break;
			case 1:
				if (THIDXn[BRCindex] <= 3)
				{
					if (hcode[hcode_index].mcode < 4)
						result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
					else
						result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC1[THIDXn[BRCindex]];
				}
				else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL1[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
				break;
			case 2:
				if (THIDXn[BRCindex] <= 5)
				{
					if (hcode[hcode_index].mcode < 6)
						result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
					else
						result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC2[THIDXn[BRCindex]];
				}
				else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL2[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
				break;
			case 3:
				if (THIDXn[BRCindex] <= 6)
				{
					if (hcode[hcode_index].mcode < 9)
						result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
					else
						result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC3[THIDXn[BRCindex]];
				}
				else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL3[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
				break;
			case 4:
				if (THIDXn[BRCindex] <= 8)
				{
					if (hcode[hcode_index].mcode < 15)
						result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
					else
						result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC4[THIDXn[BRCindex]];
				}
				else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL4[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
				break;
			default: printf("UNHANDLED CASE\n"); exit(-1);
			}
			hcode_index++;
		}
		BRCindex++;
	} while (hcode_index < NQ);
}
