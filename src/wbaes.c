#include "wbaes.h"

void printstate(unsigned char * in)
{
    for(int i = 0; i < 16; i++)
    {
        printf("%.2X", in[i]);
    }
    printf("\n");
    return ;
}

void wbaes_gen(u8 key[16])
{
    u8 expandedKey[176];
    M8 L[9][16];
    M8 L_inv[9][16];
    M8 MaskMat_R1[16];
    M8 MaskMat_R1_inv[16];
    M8 MaskMat_R8[16];
    M8 MaskMat_R8_inv[16];
    M8 MaskMat_R9[16];
    M8 MaskMat_R9_inv[16];
    M32 MB[9][4];
    M32 MB_inv[9][4];
    u32 Tyi[4][256];
    M32 Out_L[9][4];
    M32 Out_Mask_R1[4];
    M32 Out_Mask_R8[4];
    M32 Out_Mask_R9[4];
    u8 temp_u8;
    u32 temp_u32;
    u32 temp_mask;
    u8 temp_u8_left;
    u8 temp_u8_right;
    int i, j, x, k, r, y;
    static u8 nibble[16] = {0x01, 0x02, 0x0C, 0x05, 0x07, 0x08, 0x0A, 0x0F, 0x04, 0x0D, 0x0B, 0x0E, 0x09, 0x06, 0x00, 0x03};
    static u8 nibble_inv[16] = {0x0e, 0x00, 0x01, 0x0f, 0x08, 0x03, 0x0d, 0x04, 0x05, 0x0c, 0x06, 0x0a, 0x02, 0x09, 0x0b, 0x07}; 
    int columnindex[]={0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
    int shiftindex[]={0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11};
    int shiftbit[]={24, 16, 8, 0};

    expandKey (key, expandedKey);

    for(i = 0; i < 9; i++)
    {
        for(j = 0; j < 16; j++)
        {
            genMatpairM8(&L[i][j], &L_inv[i][j]);
        }
    }

    for(j = 0; j < 16; j++)
    {
        genMatpairM8(&MaskMat_R1[j], &MaskMat_R1_inv[j]);
        genMatpairM8(&MaskMat_R8[j], &MaskMat_R8_inv[j]);
        genMatpairM8(&MaskMat_R9[j], &MaskMat_R9_inv[j]);
    }
    for(i = 0; i < 9; i++)
    {
        for(j = 0; j < 4; j++)
        {
            genMatpairM32(&MB[i][j], &MB_inv[i][j]);
        }
    }

    for (x = 0; x < 256; x++)
    {
      Tyi[0][x] = (gMul(2, x) << 24) | (x << 16) | (x << 8) | gMul(3, x);
      Tyi[1][x] = (gMul(3, x) << 24) | (gMul(2, x) << 16) | (x << 8) | x;
      Tyi[2][x] = (x << 24) | (gMul(3, x) << 16) | (gMul(2, x) << 8) | x;
      Tyi[3][x] = (x << 24) | (x << 16) | (gMul(3, x) << 8) | gMul(2, x);
    }

    for(i = 0; i < 9; i++)
    {
        for(j = 0; j < 4; j++)
        {
            MatrixcomM8to32(L[i][4 * j], L[i][4 * j + 1], L[i][4 * j + 2], L[i][4 * j + 3], &Out_L[i][j]);
        }
    }
    
    for(j = 0; j < 4; j++)
    {
        MatrixcomM8to32(MaskMat_R1[4 * j], MaskMat_R1[4 * j + 1], MaskMat_R1[4 * j + 2], MaskMat_R1[4 * j + 3], &Out_Mask_R1[j]);
        MatrixcomM8to32(MaskMat_R8[4 * j], MaskMat_R8[4 * j + 1], MaskMat_R8[4 * j + 2], MaskMat_R8[4 * j + 3], &Out_Mask_R8[j]);
        MatrixcomM8to32(MaskMat_R9[4 * j], MaskMat_R9[4 * j + 1], MaskMat_R9[4 * j + 2], MaskMat_R9[4 * j + 3], &Out_Mask_R9[j]);
    }
    
    //Round 1
    shiftRows (expandedKey + 16 * 0);
    InitRandom(((unsigned int)time(NULL)));
    for(j = 0; j < 16; j++)//type_II
    {
        for(x = 0; x < 256; x++)
        {
            temp_u8 = SBox[x ^ expandedKey[16 * 0 + j]];
            temp_u32 = Tyi[j % 4][temp_u8];
            temp_mask = cus_random();
            temp_u32 ^= temp_mask;
            temp_u32 = MatMulNumM32(MB[0][columnindex[j]], temp_u32);
            TypeII_MO_R1[j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        
            //mask table
            temp_u32 = MatMulNumM32(Out_Mask_R1[columnindex[j]], temp_mask);
            TypeII_MO_R1_Mask[j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }
    for(j = 0; j < 16; j++)//type_III
    {
        for(x = 0; x < 256; x++)
        {
            temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]); 
            temp_u32 = temp_u8;
            temp_u32 = temp_u32 << shiftbit[j % 4];
            temp_u32 = MatMulNumM32(MB_inv[0][columnindex[j]], temp_u32);
            temp_u32 = MatMulNumM32(Out_L[0][columnindex[j]], temp_u32);
            TypeIII[0][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }

    //Round 2
    shiftRows (expandedKey + 16 * 1);
    for(j = 0; j < 16; j++)
    {
        for(x = 0; x < 256; x++)
        {
            for(y = 0; y < 256; y++)
            {
                temp_u8_left = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
                temp_u8_right = (nibble_inv[(y & 0xf0) >> 4] << 4) | (nibble_inv[(y & 0x0f)]);
                
                temp_u8_left = MatMulNumM8(L_inv[0][shiftindex[j]], temp_u8_left);
                temp_u8_right = MatMulNumM8(MaskMat_R1_inv[shiftindex[j]], temp_u8_right);
                temp_u8 = temp_u8_left ^ temp_u8_right;
                temp_u8 = SBox[temp_u8 ^ expandedKey[16 * 1 + j]];
                temp_u32 = Tyi[j % 4][temp_u8];
                temp_mask = cus_random();
                temp_u32 ^= temp_mask;
                temp_u32 = MatMulNumM32(MB[1][columnindex[j]], temp_u32);
                TypeII_MIMO_R2[j][x][y] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);

                //mask table
                temp_u32 = MatMulNumM32(MB[1][columnindex[j]], temp_mask);
                TypeII_MIMO_R2_Mask[j][x][y] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    }
    for(j = 0; j < 16; j++)//type_III
    {
        for(x = 0; x < 256; x++)
        {
            temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
            temp_u32 = temp_u8;
            temp_u32 = temp_u32 << shiftbit[j % 4];
            temp_u32 = MatMulNumM32(MB_inv[1][columnindex[j]], temp_u32);
            temp_u32 = MatMulNumM32(Out_L[1][columnindex[j]], temp_u32);
            TypeIII[1][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }

    //Round 3 - 7
    for (i = 2; i < 7; i++)//Type_II
    {
        shiftRows (expandedKey + 16 * i);
        for(j = 0; j < 16; j++)
        {
            for(x = 0; x < 256; x++)
            {
                temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
                temp_u8 = MatMulNumM8(L_inv[i - 1][shiftindex[j]], temp_u8);
                temp_u8 = SBox[temp_u8 ^ expandedKey[16 * i + j]];
                temp_u32 = Tyi[j % 4][temp_u8];
                temp_u32 = MatMulNumM32(MB[i][columnindex[j]], temp_u32);
                TypeII[i - 2][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    
        for(j = 0; j < 16; j++)//type_III
        {
            for(x = 0; x < 256; x++)
            {
                temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
                temp_u32 = temp_u8;
                temp_u32 = temp_u32 << shiftbit[j % 4];
                temp_u32 = MatMulNumM32(MB_inv[i][columnindex[j]], temp_u32);
                temp_u32 = MatMulNumM32(Out_L[i][columnindex[j]], temp_u32);
                TypeIII[i][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    }

    //Round 8
    shiftRows (expandedKey + 16 * 7);
    for(j = 0; j < 16; j++)
    {
        for(x = 0; x < 256; x++)
        {
            temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
            temp_u8 = MatMulNumM8(L_inv[6][shiftindex[j]], temp_u8);
            temp_u8 = SBox[temp_u8 ^ expandedKey[16 * 7 + j]];
            temp_u32 = Tyi[j % 4][temp_u8];
            temp_mask = cus_random();
            temp_u32 ^= temp_mask;
            temp_u32 = MatMulNumM32(MB[7][columnindex[j]], temp_u32);
            TypeII_MO_R8[j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        
            //mask table
            temp_u32 = MatMulNumM32(Out_Mask_R8[columnindex[j]], temp_mask);
            TypeII_MO_R8_Mask[j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }
    for(j = 0; j < 16; j++)//type_III
    {
        for(x = 0; x < 256; x++)
        {
            temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]); 
            temp_u32 = temp_u8;
            temp_u32 = temp_u32 << shiftbit[j % 4];
            temp_u32 = MatMulNumM32(MB_inv[7][columnindex[j]], temp_u32);
            temp_u32 = MatMulNumM32(Out_L[7][columnindex[j]], temp_u32);
            TypeIII[7][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }

    //Round 9
    shiftRows (expandedKey + 16 * 8);
    for(j = 0; j < 16; j++)
    {
        for(x = 0; x < 256; x++)
        {
            for(y = 0; y < 256; y++)
            {
                temp_u8_left = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
                temp_u8_right = (nibble_inv[(y & 0xf0) >> 4] << 4) | (nibble_inv[(y & 0x0f)]);
                
                temp_u8_left = MatMulNumM8(L_inv[7][shiftindex[j]], temp_u8_left);
                temp_u8_right = MatMulNumM8(MaskMat_R8_inv[shiftindex[j]], temp_u8_right);
                temp_u8 = temp_u8_left ^ temp_u8_right;
                temp_u8 = SBox[temp_u8 ^ expandedKey[16 * 8 + j]];
                temp_u32 = Tyi[j % 4][temp_u8];
                temp_mask = cus_random();
                temp_u32 ^= temp_mask;
                temp_u32 = MatMulNumM32(MB[8][columnindex[j]], temp_u32);
                TypeII_MIMO_R9[j][x][y] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);

                //mask table
                temp_u32 = MatMulNumM32(Out_Mask_R9[columnindex[j]], temp_mask);
                TypeII_MIMO_R9_Mask[j][x][y] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    }
    for(j = 0; j < 16; j++)//type_III
    {
        for(x = 0; x < 256; x++)
        {
            temp_u8 = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
            temp_u32 = temp_u8;
            temp_u32 = temp_u32 << shiftbit[j % 4];
            temp_u32 = MatMulNumM32(MB_inv[8][columnindex[j]], temp_u32);
            temp_u32 = MatMulNumM32(Out_L[8][columnindex[j]], temp_u32);
            TypeIII[8][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }

    //Round 10
    shiftRows (expandedKey + 16 * 9);
    for(j = 0; j < 16; j++)//type_II
    {
        for(x = 0; x < 256; x++)
        {
            for(y = 0; y < 256; y++)
            {
                temp_u8_left = (nibble_inv[(x & 0xf0) >> 4] << 4) | (nibble_inv[(x & 0x0f)]);
                temp_u8_right = (nibble_inv[(y & 0xf0) >> 4] << 4) | (nibble_inv[(y & 0x0f)]);

                temp_u8_left = MatMulNumM8(L_inv[8][shiftindex[j]], temp_u8_left);
                temp_u8_right = MatMulNumM8(MaskMat_R9_inv[shiftindex[j]], temp_u8_right);
                temp_u8 = temp_u8_left ^ temp_u8_right;
                temp_u8 = SBox[temp_u8 ^ expandedKey[16 * 9 + j]];
                TypeV_MI[j][x][y] = temp_u8 ^ expandedKey[16 * 10 + j];
            }
        }
    }

    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 4; j++)
        {
            for(k = 0; k < 3; k++)
            {
                for(r = 0; r < 8; r++)
                {
                    for (x = 0; x < 16; x++)
                    {
                        for (y = 0; y < 16; y++)
                        {
                            TypeIV_II[i][j][k][r][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                            TypeIV_III[i][j][k][r][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                        }
                    }
                }
            }
        }
    }
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 8; k++)
            {
                for (x = 0; x < 16; x++)
                {
                    for (y = 0; y < 16; y++)
                    {
                        TypeIV_IIM_R1[i][j][k][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                        TypeIV_IIM_R8[i][j][k][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                        TypeIV_IIM_R9[i][j][k][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                    }
                }
            }
        }
    }
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            for (k = 0; k < 8; k++)
            {
                for (x = 0; x < 16; x++)
                {
                    for (y = 0; y < 16; y++)
                    {
                        TypeIV_IIM_R2[i][j][k][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                    }
                }
            }
        }
    }
}

void wbaes_encrypt(u8 input[16], u8 output[16])
{
    u32 a, b, c, d, ab, cd, aa, bb, cc, dd;
    u8 state[16];
    u8 mask[16];
    int i, j;

    for(i = 0; i < 16; i++)
    {
        state[i] = input[i];
    }
    // Round 1
    shiftRows (state);
    for (j = 0; j < 4; j++)
    {
        a = TypeII_MO_R1[4*j + 0][state[4*j + 0]];
        b = TypeII_MO_R1[4*j + 1][state[4*j + 1]];
        c = TypeII_MO_R1[4*j + 2][state[4*j + 2]];
        d = TypeII_MO_R1[4*j + 3][state[4*j + 3]];

        aa = TypeII_MO_R1_Mask[4*j + 0][state[4*j + 0]];
        bb = TypeII_MO_R1_Mask[4*j + 1][state[4*j + 1]];
        cc = TypeII_MO_R1_Mask[4*j + 2][state[4*j + 2]];
        dd = TypeII_MO_R1_Mask[4*j + 3][state[4*j + 3]];
        // 
        ab = (TypeIV_II[0][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_II[0][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_II[0][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) | (TypeIV_II[0][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_II[0][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_II[0][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_II[0][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_II[0][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_II[0][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_II[0][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_II[0][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) | (TypeIV_II[0][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_II[0][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_II[0][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_II[0][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_II[0][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_II[0][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_II[0][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_II[0][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_II[0][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_II[0][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_II[0][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_II[0][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_II[0][j][2][7][ab & 0xf][cd & 0xf];
        // 
        ab = (TypeIV_IIM_R1[j][0][0][(aa >> 28) & 0xf][(bb >> 28) & 0xf] << 28) | (TypeIV_IIM_R1[j][0][1][(aa >> 24) & 0xf][(bb >> 24) & 0xf] << 24) | (TypeIV_IIM_R1[j][0][2][(aa >> 20) & 0xf][(bb >> 20) & 0xf] << 20) | (TypeIV_IIM_R1[j][0][3][(aa >> 16) & 0xf][(bb >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R1[j][0][4][(aa >> 12) & 0xf][(bb >> 12) & 0xf] << 12) | (TypeIV_IIM_R1[j][0][5][(aa >> 8) & 0xf][(bb >> 8) & 0xf] << 8) | (TypeIV_IIM_R1[j][0][6][(aa >> 4) & 0xf][(bb >> 4) & 0xf] << 4) | TypeIV_IIM_R1[j][0][7][aa & 0xf][bb & 0xf];
        
        cd = (TypeIV_IIM_R1[j][1][0][(cc >> 28) & 0xf][(dd >> 28) & 0xf] << 28) | (TypeIV_IIM_R1[j][1][1][(cc >> 24) & 0xf][(dd >> 24) & 0xf] << 24) | (TypeIV_IIM_R1[j][1][2][(cc >> 20) & 0xf][(dd >> 20) & 0xf] << 20) | (TypeIV_IIM_R1[j][1][3][(cc >> 16) & 0xf][(dd >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R1[j][1][4][(cc >> 12) & 0xf][(dd >> 12) & 0xf] << 12) | (TypeIV_IIM_R1[j][1][5][(cc >> 8) & 0xf][(dd >> 8) & 0xf] << 8) | (TypeIV_IIM_R1[j][1][6][(cc >> 4) & 0xf][(dd >> 4) & 0xf] << 4) | TypeIV_IIM_R1[j][1][7][cc & 0xf][dd & 0xf];
        
        mask[4*j + 0] = (TypeIV_IIM_R1[j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_IIM_R1[j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        mask[4*j + 1] = (TypeIV_IIM_R1[j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_IIM_R1[j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        mask[4*j + 2] = (TypeIV_IIM_R1[j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_IIM_R1[j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        mask[4*j + 3] = (TypeIV_IIM_R1[j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_IIM_R1[j][2][7][ab & 0xf][cd & 0xf];
        // 
        a = TypeIII[0][4*j + 0][state[4*j + 0]];
        b = TypeIII[0][4*j + 1][state[4*j + 1]];
        c = TypeIII[0][4*j + 2][state[4*j + 2]];
        d = TypeIII[0][4*j + 3][state[4*j + 3]];

        ab = (TypeIV_III[0][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_III[0][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_III[0][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) | (TypeIV_III[0][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_III[0][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_III[0][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_III[0][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_III[0][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_III[0][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_III[0][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_III[0][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) | (TypeIV_III[0][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_III[0][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_III[0][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_III[0][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_III[0][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_III[0][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_III[0][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_III[0][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_III[0][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_III[0][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_III[0][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_III[0][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_III[0][j][2][7][ab & 0xf][cd & 0xf];
    }

    // Round 2
    shiftRows (state);
    shiftRows (mask);
    for (j = 0; j < 4; j++)
    {
        a = TypeII_MIMO_R2[4*j + 0][state[4*j + 0]][mask[4*j + 0]];
        b = TypeII_MIMO_R2[4*j + 1][state[4*j + 1]][mask[4*j + 1]];
        c = TypeII_MIMO_R2[4*j + 2][state[4*j + 2]][mask[4*j + 2]];
        d = TypeII_MIMO_R2[4*j + 3][state[4*j + 3]][mask[4*j + 3]];

        aa = TypeII_MIMO_R2_Mask[4*j + 0][state[4*j + 0]][mask[4*j + 0]];
        bb = TypeII_MIMO_R2_Mask[4*j + 1][state[4*j + 1]][mask[4*j + 1]];
        cc = TypeII_MIMO_R2_Mask[4*j + 2][state[4*j + 2]][mask[4*j + 2]];
        dd = TypeII_MIMO_R2_Mask[4*j + 3][state[4*j + 3]][mask[4*j + 3]];
        // 
        ab = (TypeIV_II[1][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_II[1][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_II[1][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) | (TypeIV_II[1][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_II[1][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_II[1][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_II[1][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_II[1][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_II[1][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_II[1][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_II[1][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) | (TypeIV_II[1][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_II[1][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_II[1][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_II[1][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_II[1][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_II[1][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_II[1][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_II[1][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_II[1][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_II[1][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_II[1][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_II[1][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_II[1][j][2][7][ab & 0xf][cd & 0xf];
        // 
        ab = (TypeIV_IIM_R2[j][0][0][(aa >> 28) & 0xf][(bb >> 28) & 0xf] << 28) | (TypeIV_IIM_R2[j][0][1][(aa >> 24) & 0xf][(bb >> 24) & 0xf] << 24) | (TypeIV_IIM_R2[j][0][2][(aa >> 20) & 0xf][(bb >> 20) & 0xf] << 20) | (TypeIV_IIM_R2[j][0][3][(aa >> 16) & 0xf][(bb >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R2[j][0][4][(aa >> 12) & 0xf][(bb >> 12) & 0xf] << 12) | (TypeIV_IIM_R2[j][0][5][(aa >> 8) & 0xf][(bb >> 8) & 0xf] << 8) | (TypeIV_IIM_R2[j][0][6][(aa >> 4) & 0xf][(bb >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][0][7][aa & 0xf][bb & 0xf];
        
        cd = (TypeIV_IIM_R2[j][1][0][(cc >> 28) & 0xf][(dd >> 28) & 0xf] << 28) | (TypeIV_IIM_R2[j][1][1][(cc >> 24) & 0xf][(dd >> 24) & 0xf] << 24) | (TypeIV_IIM_R2[j][1][2][(cc >> 20) & 0xf][(dd >> 20) & 0xf] << 20) | (TypeIV_IIM_R2[j][1][3][(cc >> 16) & 0xf][(dd >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R2[j][1][4][(cc >> 12) & 0xf][(dd >> 12) & 0xf] << 12) | (TypeIV_IIM_R2[j][1][5][(cc >> 8) & 0xf][(dd >> 8) & 0xf] << 8) | (TypeIV_IIM_R2[j][1][6][(cc >> 4) & 0xf][(dd >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][1][7][cc & 0xf][dd & 0xf];
        
        mask[4*j + 0] = (TypeIV_IIM_R2[j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_IIM_R2[j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        mask[4*j + 1] = (TypeIV_IIM_R2[j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_IIM_R2[j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        mask[4*j + 2] = (TypeIV_IIM_R2[j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_IIM_R2[j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        mask[4*j + 3] = (TypeIV_IIM_R2[j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][2][7][ab & 0xf][cd & 0xf];

        state[4*j + 0] = (TypeIV_IIM_R2[j][3][0][(state[4*j + 0] >> 4) & 0xf][(mask[4*j + 0] >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][3][1][state[4*j + 0] & 0xf][mask[4*j + 0] & 0xf];
        state[4*j + 1] = (TypeIV_IIM_R2[j][3][2][(state[4*j + 1] >> 4) & 0xf][(mask[4*j + 1] >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][3][3][state[4*j + 1] & 0xf][mask[4*j + 1] & 0xf];
        state[4*j + 2] = (TypeIV_IIM_R2[j][3][4][(state[4*j + 2] >> 4) & 0xf][(mask[4*j + 2] >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][3][5][state[4*j + 2] & 0xf][mask[4*j + 2] & 0xf];
        state[4*j + 3] = (TypeIV_IIM_R2[j][3][6][(state[4*j + 3] >> 4) & 0xf][(mask[4*j + 3] >> 4) & 0xf] << 4) | TypeIV_IIM_R2[j][3][7][state[4*j + 3] & 0xf][mask[4*j + 3] & 0xf];

        a = TypeIII[1][4*j + 0][state[4*j + 0]];
        b = TypeIII[1][4*j + 1][state[4*j + 1]];
        c = TypeIII[1][4*j + 2][state[4*j + 2]];
        d = TypeIII[1][4*j + 3][state[4*j + 3]];

        ab = (TypeIV_III[1][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_III[1][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_III[1][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) |(TypeIV_III[1][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_III[1][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_III[1][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_III[1][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_III[1][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_III[1][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_III[1][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_III[1][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) |(TypeIV_III[1][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_III[1][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_III[1][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_III[1][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_III[1][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_III[1][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_III[1][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_III[1][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_III[1][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_III[1][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_III[1][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_III[1][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_III[1][j][2][7][ab & 0xf][cd & 0xf];
    }

    // Round 3 - 7
    for (i = 2; i < 7; i++)
    {
        shiftRows (state);
        for (j = 0; j < 4; j++)
        {
            a = TypeII[i - 2][4*j + 0][state[4*j + 0]];
            b = TypeII[i - 2][4*j + 1][state[4*j + 1]];
            c = TypeII[i - 2][4*j + 2][state[4*j + 2]];
            d = TypeII[i - 2][4*j + 3][state[4*j + 3]];

            ab = (TypeIV_II[i][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_II[i][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_II[i][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) |(TypeIV_II[i][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
            (TypeIV_II[i][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_II[i][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_II[i][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_II[i][j][0][7][a & 0xf][b & 0xf];
            
            cd = (TypeIV_II[i][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_II[i][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_II[i][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) |(TypeIV_II[i][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
            (TypeIV_II[i][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_II[i][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_II[i][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_II[i][j][1][7][c & 0xf][d & 0xf];
            
            state[4*j + 0] = (TypeIV_II[i][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_II[i][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
            state[4*j + 1] = (TypeIV_II[i][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_II[i][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
            state[4*j + 2] = (TypeIV_II[i][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_II[i][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
            state[4*j + 3] = (TypeIV_II[i][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_II[i][j][2][7][ab & 0xf][cd & 0xf];

            a = TypeIII[i][4*j + 0][state[4*j + 0]];
            b = TypeIII[i][4*j + 1][state[4*j + 1]];
            c = TypeIII[i][4*j + 2][state[4*j + 2]];
            d = TypeIII[i][4*j + 3][state[4*j + 3]];

            ab = (TypeIV_III[i][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_III[i][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_III[i][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) |(TypeIV_III[i][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
            (TypeIV_III[i][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_III[i][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_III[i][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_III[i][j][0][7][a & 0xf][b & 0xf];
            
            cd = (TypeIV_III[i][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_III[i][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_III[i][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) |(TypeIV_III[i][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
            (TypeIV_III[i][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_III[i][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_III[i][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_III[i][j][1][7][c & 0xf][d & 0xf];
            
            state[4*j + 0] = (TypeIV_III[i][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_III[i][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
            state[4*j + 1] = (TypeIV_III[i][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_III[i][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
            state[4*j + 2] = (TypeIV_III[i][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_III[i][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
            state[4*j + 3] = (TypeIV_III[i][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_III[i][j][2][7][ab & 0xf][cd & 0xf];
        }
    }

    // Round 8
    shiftRows (state);
    for (j = 0; j < 4; j++)
    {
        a = TypeII_MO_R8[4*j + 0][state[4*j + 0]];
        b = TypeII_MO_R8[4*j + 1][state[4*j + 1]];
        c = TypeII_MO_R8[4*j + 2][state[4*j + 2]];
        d = TypeII_MO_R8[4*j + 3][state[4*j + 3]];

        aa = TypeII_MO_R8_Mask[4*j + 0][state[4*j + 0]];
        bb = TypeII_MO_R8_Mask[4*j + 1][state[4*j + 1]];
        cc = TypeII_MO_R8_Mask[4*j + 2][state[4*j + 2]];
        dd = TypeII_MO_R8_Mask[4*j + 3][state[4*j + 3]];
        // 
        ab = (TypeIV_II[7][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_II[7][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_II[7][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) | (TypeIV_II[7][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_II[7][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_II[7][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_II[7][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_II[7][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_II[7][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_II[7][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_II[7][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) | (TypeIV_II[7][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_II[7][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_II[7][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_II[7][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_II[7][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_II[7][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_II[7][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_II[7][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_II[7][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_II[7][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_II[7][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_II[7][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_II[7][j][2][7][ab & 0xf][cd & 0xf];
        // 
        ab = (TypeIV_IIM_R8[j][0][0][(aa >> 28) & 0xf][(bb >> 28) & 0xf] << 28) | (TypeIV_IIM_R8[j][0][1][(aa >> 24) & 0xf][(bb >> 24) & 0xf] << 24) | (TypeIV_IIM_R8[j][0][2][(aa >> 20) & 0xf][(bb >> 20) & 0xf] << 20) | (TypeIV_IIM_R8[j][0][3][(aa >> 16) & 0xf][(bb >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R8[j][0][4][(aa >> 12) & 0xf][(bb >> 12) & 0xf] << 12) | (TypeIV_IIM_R8[j][0][5][(aa >> 8) & 0xf][(bb >> 8) & 0xf] << 8) | (TypeIV_IIM_R8[j][0][6][(aa >> 4) & 0xf][(bb >> 4) & 0xf] << 4) | TypeIV_IIM_R8[j][0][7][aa & 0xf][bb & 0xf];
        
        cd = (TypeIV_IIM_R8[j][1][0][(cc >> 28) & 0xf][(dd >> 28) & 0xf] << 28) | (TypeIV_IIM_R8[j][1][1][(cc >> 24) & 0xf][(dd >> 24) & 0xf] << 24) | (TypeIV_IIM_R8[j][1][2][(cc >> 20) & 0xf][(dd >> 20) & 0xf] << 20) | (TypeIV_IIM_R8[j][1][3][(cc >> 16) & 0xf][(dd >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R8[j][1][4][(cc >> 12) & 0xf][(dd >> 12) & 0xf] << 12) | (TypeIV_IIM_R8[j][1][5][(cc >> 8) & 0xf][(dd >> 8) & 0xf] << 8) | (TypeIV_IIM_R8[j][1][6][(cc >> 4) & 0xf][(dd >> 4) & 0xf] << 4) | TypeIV_IIM_R8[j][1][7][cc & 0xf][dd & 0xf];
        
        mask[4*j + 0] = (TypeIV_IIM_R8[j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_IIM_R8[j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        mask[4*j + 1] = (TypeIV_IIM_R8[j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_IIM_R8[j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        mask[4*j + 2] = (TypeIV_IIM_R8[j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_IIM_R8[j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        mask[4*j + 3] = (TypeIV_IIM_R8[j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_IIM_R8[j][2][7][ab & 0xf][cd & 0xf];
        // 
        a = TypeIII[7][4*j + 0][state[4*j + 0]];
        b = TypeIII[7][4*j + 1][state[4*j + 1]];
        c = TypeIII[7][4*j + 2][state[4*j + 2]];
        d = TypeIII[7][4*j + 3][state[4*j + 3]];

        ab = (TypeIV_III[7][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_III[7][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_III[7][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) | (TypeIV_III[7][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_III[7][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_III[7][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_III[7][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_III[7][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_III[7][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_III[7][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_III[7][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) | (TypeIV_III[7][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_III[7][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_III[7][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_III[7][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_III[7][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_III[7][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_III[7][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_III[7][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_III[7][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_III[7][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_III[7][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_III[7][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_III[7][j][2][7][ab & 0xf][cd & 0xf];
    }

    // Round 9
    shiftRows (state);
    shiftRows (mask);
    for (j = 0; j < 4; j++)
    {
        a = TypeII_MIMO_R9[4*j + 0][state[4*j + 0]][mask[4*j + 0]];
        b = TypeII_MIMO_R9[4*j + 1][state[4*j + 1]][mask[4*j + 1]];
        c = TypeII_MIMO_R9[4*j + 2][state[4*j + 2]][mask[4*j + 2]];
        d = TypeII_MIMO_R9[4*j + 3][state[4*j + 3]][mask[4*j + 3]];

        aa = TypeII_MIMO_R9_Mask[4*j + 0][state[4*j + 0]][mask[4*j + 0]];
        bb = TypeII_MIMO_R9_Mask[4*j + 1][state[4*j + 1]][mask[4*j + 1]];
        cc = TypeII_MIMO_R9_Mask[4*j + 2][state[4*j + 2]][mask[4*j + 2]];
        dd = TypeII_MIMO_R9_Mask[4*j + 3][state[4*j + 3]][mask[4*j + 3]];
        // 
        ab = (TypeIV_II[8][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_II[8][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_II[8][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) | (TypeIV_II[8][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_II[8][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_II[8][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_II[8][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_II[8][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_II[8][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_II[8][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_II[8][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) | (TypeIV_II[8][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_II[8][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_II[8][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_II[8][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_II[8][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_II[8][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_II[8][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_II[8][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_II[8][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_II[8][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_II[8][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_II[8][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_II[8][j][2][7][ab & 0xf][cd & 0xf];
        // 
        ab = (TypeIV_IIM_R9[j][0][0][(aa >> 28) & 0xf][(bb >> 28) & 0xf] << 28) | (TypeIV_IIM_R9[j][0][1][(aa >> 24) & 0xf][(bb >> 24) & 0xf] << 24) | (TypeIV_IIM_R9[j][0][2][(aa >> 20) & 0xf][(bb >> 20) & 0xf] << 20) | (TypeIV_IIM_R9[j][0][3][(aa >> 16) & 0xf][(bb >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R9[j][0][4][(aa >> 12) & 0xf][(bb >> 12) & 0xf] << 12) | (TypeIV_IIM_R9[j][0][5][(aa >> 8) & 0xf][(bb >> 8) & 0xf] << 8) | (TypeIV_IIM_R9[j][0][6][(aa >> 4) & 0xf][(bb >> 4) & 0xf] << 4) | TypeIV_IIM_R9[j][0][7][aa & 0xf][bb & 0xf];
        
        cd = (TypeIV_IIM_R9[j][1][0][(cc >> 28) & 0xf][(dd >> 28) & 0xf] << 28) | (TypeIV_IIM_R9[j][1][1][(cc >> 24) & 0xf][(dd >> 24) & 0xf] << 24) | (TypeIV_IIM_R9[j][1][2][(cc >> 20) & 0xf][(dd >> 20) & 0xf] << 20) | (TypeIV_IIM_R9[j][1][3][(cc >> 16) & 0xf][(dd >> 16) & 0xf] << 16) |\
        (TypeIV_IIM_R9[j][1][4][(cc >> 12) & 0xf][(dd >> 12) & 0xf] << 12) | (TypeIV_IIM_R9[j][1][5][(cc >> 8) & 0xf][(dd >> 8) & 0xf] << 8) | (TypeIV_IIM_R9[j][1][6][(cc >> 4) & 0xf][(dd >> 4) & 0xf] << 4) | TypeIV_IIM_R9[j][1][7][cc & 0xf][dd & 0xf];
        
        mask[4*j + 0] = (TypeIV_IIM_R9[j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_IIM_R9[j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        mask[4*j + 1] = (TypeIV_IIM_R9[j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_IIM_R9[j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        mask[4*j + 2] = (TypeIV_IIM_R9[j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_IIM_R9[j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        mask[4*j + 3] = (TypeIV_IIM_R9[j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_IIM_R9[j][2][7][ab & 0xf][cd & 0xf];

        a = TypeIII[8][4*j + 0][state[4*j + 0]];
        b = TypeIII[8][4*j + 1][state[4*j + 1]];
        c = TypeIII[8][4*j + 2][state[4*j + 2]];
        d = TypeIII[8][4*j + 3][state[4*j + 3]];

        ab = (TypeIV_III[8][j][0][0][(a >> 28) & 0xf][(b >> 28) & 0xf] << 28) | (TypeIV_III[8][j][0][1][(a >> 24) & 0xf][(b >> 24) & 0xf] << 24) | (TypeIV_III[8][j][0][2][(a >> 20) & 0xf][(b >> 20) & 0xf] << 20) |(TypeIV_III[8][j][0][3][(a >> 16) & 0xf][(b >> 16) & 0xf] << 16) |\
        (TypeIV_III[8][j][0][4][(a >> 12) & 0xf][(b >> 12) & 0xf] << 12) | (TypeIV_III[8][j][0][5][(a >> 8) & 0xf][(b >> 8) & 0xf] << 8) | (TypeIV_III[8][j][0][6][(a >> 4) & 0xf][(b >> 4) & 0xf] << 4) | TypeIV_III[8][j][0][7][a & 0xf][b & 0xf];
        
        cd = (TypeIV_III[8][j][1][0][(c >> 28) & 0xf][(d >> 28) & 0xf] << 28) | (TypeIV_III[8][j][1][1][(c >> 24) & 0xf][(d >> 24) & 0xf] << 24) | (TypeIV_III[8][j][1][2][(c >> 20) & 0xf][(d >> 20) & 0xf] << 20) |(TypeIV_III[8][j][1][3][(c >> 16) & 0xf][(d >> 16) & 0xf] << 16) |\
        (TypeIV_III[8][j][1][4][(c >> 12) & 0xf][(d >> 12) & 0xf] << 12) | (TypeIV_III[8][j][1][5][(c >> 8) & 0xf][(d >> 8) & 0xf] << 8) | (TypeIV_III[8][j][1][6][(c >> 4) & 0xf][(d >> 4) & 0xf] << 4) | TypeIV_III[8][j][1][7][c & 0xf][d & 0xf];
        
        state[4*j + 0] = (TypeIV_III[8][j][2][0][(ab >> 28) & 0xf][(cd >> 28) & 0xf] << 4) | TypeIV_III[8][j][2][1][(ab >> 24) & 0xf][(cd >> 24) & 0xf];
        state[4*j + 1] = (TypeIV_III[8][j][2][2][(ab >> 20) & 0xf][(cd >> 20) & 0xf] << 4) | TypeIV_III[8][j][2][3][(ab >> 16) & 0xf][(cd >> 16) & 0xf];
        state[4*j + 2] = (TypeIV_III[8][j][2][4][(ab >> 12) & 0xf][(cd >> 12) & 0xf] << 4) | TypeIV_III[8][j][2][5][(ab >> 8) & 0xf][(cd >> 8) & 0xf];
        state[4*j + 3] = (TypeIV_III[8][j][2][6][(ab >> 4) & 0xf][(cd >> 4) & 0xf] << 4) | TypeIV_III[8][j][2][7][ab & 0xf][cd & 0xf];
    }

    // Round 10
    shiftRows (state);
    shiftRows (mask);
    for (j = 0; j < 16; j++)
    {
        output[j] = TypeV_MI[j][state[j]][mask[j]];
    }
}