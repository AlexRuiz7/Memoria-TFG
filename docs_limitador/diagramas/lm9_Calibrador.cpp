#define HwDataReadTimesPerSec 10

/**
 * @brief Lee los datos de micrófono, líneas y los niveles de hardware
 * @param count número de muestras a tomar
 * @param mic espectro de onda del micrófono
 * @param left espectro de onda de la linea izquierda
 * @param right espectro de onda de la linea derecha
 * @param lLevel nivel de hardware. Canal izquierdo
 * @param rLevel nivel de hardware. Canal derecho
 * @param status objeto de acceso a memoria compartida
 * @return (void)
 */
void readAudioData(int count, Spectrum &mic, Spectrum &left, Spectrum &right,
                   float &lLevel, float &rLevel, HardwareStatus *status)
{
    ChannelInformation m, m2, l, r;
    Spectrum mAc, lAc, rAc;

    for (int i = 0; i < count; i++)
    {
        for (int j = 0; j < HwDataReadTimesPerSec; j++)
        {
            usleep(UnSegundo / HwDataReadTimesPerSec);
            lLevel += status->inputPressure[0];
            rLevel += status->inputPressure[1];
        };
        loadSlowChannels(m, m2, l, r);
        mAc.sum(m.spectrum);
        lAc.sum(l.spectrum);
        rAc.sum(r.spectrum);
    };

    mAc.divide(count);
    lAc.divide(count);
    rAc.divide(count);

    mic = mAc;
    left = lAc;
    right = rAc;

    lLevel /= count * HwDataReadTimesPerSec;
    rLevel /= count * HwDataReadTimesPerSec;
};

/**
 * @brief Carga las calibraciones
 * @param micCal
 * @param leftCal
 * @param rightCal
 * @return (void)
 */
void loadCalibrations(Calibracion &micCal, Calibracion &leftCal, Calibracion &rightCal)
{
    micCal.leeCalibracion(0);
    leftCal.leeCalibracion(1);
    rightCal.leeCalibracion(2);
};

/**
 * @brief Guarda las calibraciones
 * @param leftCal
 * @param rightCal
 * @return (void)
 */
void saveCalibrations(Calibracion &leftCal, Calibracion &rightCal)
{
    leftCal.guardaCalibracion(1);
    rightCal.guardaCalibracion(2);
}

/**
 * @brief Programa principal
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char **argv)
{
    Spectrum mic, left, right;
    Calibracion micCal, leftCal, rightCal;
    HardwareStatus *status;
    float lLevel, rLevel;
    float micBandValue;

    ...

    loadCalibrations(micCal, leftCal, rightCal);
    int count = 10;

    ...

    printf("--------- CALIBRADOR ---------\n");
    printf("------- Leyendo Lineas -------\n");
    readAudioData(count, mic, left, right, lLevel, rLevel, status);

    float Mic1K = mic.dB[Spectrum_1Khz_BandIndex] - micCal.internalEqualization[Spectrum_1Khz_BandIndex];

    // Referencia
    leftCal.dBRef = Mic1K;
    leftCal.ref = left.energy[Spectrum_1Khz_BandIndex];

    rightCal.dBRef = Mic1K;
    rightCal.ref = right.energy[Spectrum_1Khz_BandIndex];

    // Limpiar ecualización interna
    for (int i = 0; i < Spectrum_BandCount; i++)
        leftCal.internalEqualization[i] = rightCal.internalEqualization[i] = 0;

    // Ajustamos la banda de 1K como central
    left.calibrate(leftCal);
    right.calibrate(rightCal);

    // Calculamos la ecualización interna
    for (int i = 0; i < Spectrum_BandCount; i++)
    {
        micBandValue = mic.dB[i] - micCal.internalEqualization[i];
        leftCal.internalEqualization[i] = micBandValue - left.dB[i];
        rightCal.internalEqualization[i] = micBandValue - right.dB[i];
    }

    ...

    // Finishing
    printf("--------- Linea Left ---------\n");
    leftCal.fullPrint();
    printf("--------- Linea Right --------\n");
    rightCal.fullPrint();

    saveCalibrations(leftCal, rightCal);
    return 0;
};