#ifndef GENTABLES_H
#define GENTABLES_H

#include "aes.h"
#include "WBMatrix/WBMatrix.h"

u32 TypeII_MO_R1[16][256];
u32 TypeII_MO_R1_Mask[16][256];
u32 TypeII_MO_R8[16][256];
u32 TypeII_MO_R8_Mask[16][256];
u32 TypeII_MIMO_R2[16][256][256];
u32 TypeII_MIMO_R2_Mask[16][256][256];
u32 TypeII_MIMO_R9[16][256][256];
u32 TypeII_MIMO_R9_Mask[16][256][256];
u8 TypeV_MI[16][256][256];
u32 TypeII[5][16][256];//R3 - R7 Type II
u32 TypeIII[9][16][256];//R1 - R9 Type III
u8 TypeIV_II[9][4][3][8][16][16]; 
u8 TypeIV_III[9][4][3][8][16][16];
u8 TypeIV_IIM_R1[4][3][8][16][16];
u8 TypeIV_IIM_R8[4][3][8][16][16];
u8 TypeIV_IIM_R9[4][3][8][16][16];
u8 TypeIV_IIM_R2[4][4][8][16][16];

void wbaes_gen(u8 key[16]);
void wbaes_encrypt(u8 input[16], u8 output[16]);
void printstate(unsigned char * in);

#endif // GENTABLES_H