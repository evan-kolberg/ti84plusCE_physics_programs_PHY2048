/*
 * Projectile Motion Solver for TI-84 Plus CE
 *
 * A kinematic equation solver for 2D projectile motion problems.
 * Supports solving for position, velocity, acceleration, displacement, and time
 * in both X and Y components, with automatic equation selection and visualization.
 *
 * Controls:
 *   Arrow Keys  - Navigate between cells
 *   Enter       - Begin editing selected cell
 *   0-9, ., (-) - Enter numeric values
 *   Del         - Clear selected cell (or backspace in edit mode)
 *   Clear       - Cancel edit / Exit program
 *   Mode        - Reset all values
 *   Graph       - Show full-screen trajectory plot
 */

#include <tice.h>
#include <graphx.h>
#include <keypadc.h>
#include <ti/real.h>
#include <sys/power.h>
#include <sys/timers.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define GRAVITY 9.81f
#define INACTIVITY_TIMEOUT (5 * 60)

clock_t lastActivityTime;
#define PI 3.14159265f
#define DEG_TO_RAD (PI / 180.0f)
#define RAD_TO_DEG (180.0f / PI)

#define VAR_COUNT 7
#define ROWS 7
#define COLS 2

float xVals[VAR_COUNT];
float yVals[VAR_COUNT];
bool xKnown[VAR_COUNT];
bool yKnown[VAR_COUNT];
bool xUserSet[VAR_COUNT];
bool yUserSet[VAR_COUNT];

float launchSpeed = 0;
float launchAngle = 0;
float finalSpeed = 0;
bool speedKnown = false;
bool angleKnown = false;
bool finalSpeedKnown = false;
bool speedUserSet = false;
bool angleUserSet = false;
bool finalSpeedUserSet = false;

int curRow = 0;
int curCol = 0;
char inputBuf[12];
int inputLen = 0;
bool inputMode = false;
bool hasDecimal = false;
bool isNegative = false;

const char* rowLabels[VAR_COUNT] = {"p0", "pf", "v0", "vf", "a", "d", "t"};

char xEqUsed[64] = "";
char yEqUsed[64] = "";

void floatToStr(float val, char* out) {
    real_t r = os_FloatToReal(val);
    os_RealToStr(out, &r, 8, 1, 3);
    for (int i = 0; out[i]; i++) {
        if (out[i] == '\032') out[i] = '-';
    }
}

void initData() {
    for (int i = 0; i < VAR_COUNT; i++) {
        xVals[i] = 0;
        yVals[i] = 0;
        xKnown[i] = false;
        yKnown[i] = false;
        xUserSet[i] = false;
        yUserSet[i] = false;
    }
    xVals[0] = 0;
    xKnown[0] = true;
    xUserSet[0] = true;
    xVals[4] = 0;
    xKnown[4] = true;
    xUserSet[4] = true;
    yVals[4] = -GRAVITY;
    yKnown[4] = true;
    yUserSet[4] = true;
    
    launchSpeed = 0;
    launchAngle = 0;
    finalSpeed = 0;
    speedKnown = false;
    angleKnown = false;
    finalSpeedKnown = false;
    speedUserSet = false;
    angleUserSet = false;
    finalSpeedUserSet = false;
}

void addEq(char* eqUsed, const char* eq) {
    if (strstr(eqUsed, eq) != NULL) return;
    if (eqUsed[0] == '\0') {
        strcpy(eqUsed, eq);
    } else if (strlen(eqUsed) + strlen(eq) + 2 < 60) {
        strcat(eqUsed, "|");
        strcat(eqUsed, eq);
    }
}

void trySolve(float* vals, bool* known, bool* userSet, char* eqUsed) {
    float p0 = vals[0], pf = vals[1], v0 = vals[2], vf = vals[3], a = vals[4], d = vals[5], t = vals[6];
    bool kp0 = known[0], kpf = known[1], kv0 = known[2], kvf = known[3], ka = known[4], kd = known[5], kt = known[6];
    
    for (int iter = 0; iter < 20; iter++) {
        if (!kd && kp0 && kpf) {
            d = pf - p0;
            kd = true;
            addEq(eqUsed, "d = pf - p0");
        }
        if (!kpf && kp0 && kd) {
            pf = p0 + d;
            kpf = true;
            addEq(eqUsed, "pf = p0 + d");
        }
        if (!kp0 && kpf && kd) {
            p0 = pf - d;
            kp0 = true;
            addEq(eqUsed, "p0 = pf - d");
        }
        if (!kt && kv0 && kvf && ka && a != 0) {
            t = (vf - v0) / a;
            if (t >= 0) { kt = true; addEq(eqUsed, "t = (vf - v0) / a"); }
        }
        if (!kvf && kv0 && ka && a == 0) {
            vf = v0;
            kvf = true;
            addEq(eqUsed, "vf = v0 (a=0)");
        }
        if (!kv0 && kvf && ka && a == 0) {
            v0 = vf;
            kv0 = true;
            addEq(eqUsed, "v0 = vf (a=0)");
        }
        if (!kvf && kv0 && ka && kt) {
            vf = v0 + a * t;
            kvf = true;
            addEq(eqUsed, "vf = v0 + a*t");
        }
        if (!kv0 && kvf && ka && kt) {
            v0 = vf - a * t;
            kv0 = true;
            addEq(eqUsed, "v0 = vf - a*t");
        }
        if (!ka && kv0 && kvf && kt && t != 0) {
            a = (vf - v0) / t;
            ka = true;
            addEq(eqUsed, "a = (vf - v0) / t");
        }
        if (!kd && kv0 && kt && ka) {
            d = v0 * t + 0.5f * a * t * t;
            kd = true;
            addEq(eqUsed, "d = v0*t + .5*a*t^2");
        }
        if (!kv0 && kd && kt && ka && t != 0) {
            v0 = (d - 0.5f * a * t * t) / t;
            kv0 = true;
            addEq(eqUsed, "v0 = (d - .5*a*t^2) / t");
        }
        if (!ka && kv0 && kd && kt && t != 0) {
            a = 2 * (d - v0 * t) / (t * t);
            ka = true;
            addEq(eqUsed, "a = 2(d - v0*t) / t^2");
        }
        if (!kd && kv0 && kvf && kt) {
            d = (v0 + vf) * 0.5f * t;
            kd = true;
            addEq(eqUsed, "d = (v0 + vf) * t / 2");
        }
        if (!kt && kv0 && kvf && kd && (v0 + vf) != 0) {
            t = 2 * d / (v0 + vf);
            if (t >= 0) { kt = true; addEq(eqUsed, "t = 2*d / (v0 + vf)"); }
        }
        if (!kv0 && kvf && kd && kt && t != 0) {
            v0 = 2 * d / t - vf;
            kv0 = true;
            addEq(eqUsed, "v0 = 2*d / t - vf");
        }
        if (!kvf && kv0 && kd && kt && t != 0) {
            vf = 2 * d / t - v0;
            kvf = true;
            addEq(eqUsed, "vf = 2*d / t - v0");
        }

        if (!kd && kvf && kt && ka) {
            d = vf * t - 0.5f * a * t * t;
            kd = true;
            addEq(eqUsed, "d = vf*t - .5*a*t^2");
        }
        if (!kvf && kd && kt && ka && t != 0) {
            vf = (d + 0.5f * a * t * t) / t;
            kvf = true;
            addEq(eqUsed, "vf = (d + .5*a*t^2) / t");
        }
        if (!ka && kvf && kd && kt && t != 0) {
            a = 2 * (vf * t - d) / (t * t);
            ka = true;
            addEq(eqUsed, "a = 2(vf*t - d) / t^2");
        }
        if (!kt && kvf && kd && ka) {
            if (a != 0) {
                float disc = vf * vf - 2 * a * d;
                if (disc >= 0) {
                    float t1 = (vf - sqrtf(disc)) / a;
                    float t2 = (vf + sqrtf(disc)) / a;
                    float tMax = (t1 > t2) ? t1 : t2;
                    float tMin = (t1 < t2) ? t1 : t2;
                    if (tMax > 0.001f) {
                        t = tMax;
                        kt = true;
                    } else if (tMin >= 0) {
                        t = tMin;
                        kt = true;
                    }
                    if (kt) addEq(eqUsed, "d = vf*t - .5*a*t^2");
                }
            } else if (vf != 0) {
                t = d / vf;
                if (t >= 0) { kt = true; addEq(eqUsed, "t = d / vf"); }
            }
        }
        if (!kt && kv0 && kd && ka) {
            if (a != 0) {
                float disc = v0 * v0 + 2 * a * d;
                if (disc >= 0) {
                    float t1 = (-v0 + sqrtf(disc)) / a;
                    float t2 = (-v0 - sqrtf(disc)) / a;
                    float tMax = (t1 > t2) ? t1 : t2;
                    float tMin = (t1 < t2) ? t1 : t2;
                    if (tMax > 0.001f) {
                        t = tMax;
                        kt = true;
                    } else if (tMin >= 0) {
                        t = tMin;
                        kt = true;
                    }
                    if (kt) addEq(eqUsed, "d = v0*t + .5*a*t^2");
                }
            } else if (v0 != 0) {
                t = d / v0;
                if (t >= 0) { kt = true; addEq(eqUsed, "t = d / v0"); }
            }
        }
        if (!kvf && kv0 && ka && kd) {
            float disc = v0 * v0 + 2 * a * d;
            if (disc >= 0) {
                float vfMag = sqrtf(disc);
                if (v0 != 0) vf = (v0 > 0) ? vfMag : -vfMag;
                else vf = (a * d >= 0) ? vfMag : -vfMag;
                kvf = true;
                addEq(eqUsed, "vf^2 = v0^2 + 2*a*d");
            }
        }
        if (!kv0 && kvf && ka && kd) {
            float disc = vf * vf - 2 * a * d;
            if (disc >= 0) {
                float v0Mag = sqrtf(disc);
                if (vf != 0) v0 = (vf > 0) ? v0Mag : -v0Mag;
                else v0 = (a * d <= 0) ? v0Mag : -v0Mag;
                kv0 = true;
                addEq(eqUsed, "v0^2 = vf^2 - 2*a*d");
            }
        }
        if (!ka && kv0 && kvf && kd && d != 0) {
            a = (vf * vf - v0 * v0) / (2 * d);
            ka = true;
            addEq(eqUsed, "a = (vf^2 - v0^2) / 2*d");
        }
        if (!kd && kv0 && kvf && ka && a != 0) {
            d = (vf * vf - v0 * v0) / (2 * a);
            kd = true;
            addEq(eqUsed, "d = (vf^2 - v0^2) / 2*a");
        }
    }
    
    vals[0] = p0; vals[1] = pf; vals[2] = v0; vals[3] = vf; vals[4] = a; vals[5] = d; vals[6] = t;
    
    for (int i = 0; i < VAR_COUNT; i++) {
        bool nowKnown = false;
        switch(i) {
            case 0: nowKnown = kp0; break;
            case 1: nowKnown = kpf; break;
            case 2: nowKnown = kv0; break;
            case 3: nowKnown = kvf; break;
            case 4: nowKnown = ka; break;
            case 5: nowKnown = kd; break;
            case 6: nowKnown = kt; break;
        }
        if (!userSet[i] && nowKnown) known[i] = true;
    }
}

void autoSolve() {
    xEqUsed[0] = '\0';
    yEqUsed[0] = '\0';
    
    for (int i = 0; i < VAR_COUNT; i++) {
        if (!xUserSet[i]) xKnown[i] = false;
        if (!yUserSet[i]) yKnown[i] = false;
    }
    if (!speedUserSet) speedKnown = false;
    if (!angleUserSet) angleKnown = false;
    if (!finalSpeedUserSet) finalSpeedKnown = false;
    
    if (speedUserSet && angleUserSet) {
        xVals[2] = launchSpeed * cosf(launchAngle * DEG_TO_RAD);
        yVals[2] = launchSpeed * sinf(launchAngle * DEG_TO_RAD);
        xKnown[2] = true;
        yKnown[2] = true;
    }
    
    if (finalSpeedUserSet && xKnown[3] && !yUserSet[3]) {
        float vfx = xVals[3];
        float vfy_sq = finalSpeed * finalSpeed - vfx * vfx;
        if (vfy_sq >= 0) {
            yVals[3] = -sqrtf(vfy_sq);
            yKnown[3] = true;
        }
    }
    if (finalSpeedUserSet && yKnown[3] && !xUserSet[3]) {
        float vfy = yVals[3];
        float vfx_sq = finalSpeed * finalSpeed - vfy * vfy;
        if (vfx_sq >= 0) {
            xVals[3] = sqrtf(vfx_sq);
            xKnown[3] = true;
        }
    }
    
    if (xUserSet[6] && !yUserSet[6]) {
        yVals[6] = xVals[6];
        yKnown[6] = true;
    } else if (yUserSet[6] && !xUserSet[6]) {
        xVals[6] = yVals[6];
        xKnown[6] = true;
    }
    
    for (int pass = 0; pass < 3; pass++) {
        trySolve(xVals, xKnown, xUserSet, xEqUsed);
        
        if (xKnown[6] && !yKnown[6]) {
            yVals[6] = xVals[6];
            yKnown[6] = true;
        }
        
        trySolve(yVals, yKnown, yUserSet, yEqUsed);
        
        if (yKnown[6] && !xKnown[6]) {
            xVals[6] = yVals[6];
            xKnown[6] = true;
        }
    }
    
    if (!speedKnown && xKnown[2] && yKnown[2]) {
        launchSpeed = sqrtf(xVals[2] * xVals[2] + yVals[2] * yVals[2]);
        speedKnown = true;
    }
    if (!angleKnown && xKnown[2] && yKnown[2]) {
        launchAngle = atan2f(yVals[2], xVals[2]) * RAD_TO_DEG;
        angleKnown = true;
    }
    if (!finalSpeedKnown && xKnown[3] && yKnown[3]) {
        finalSpeed = sqrtf(xVals[3] * xVals[3] + yVals[3] * yVals[3]);
        finalSpeedKnown = true;
    }
}

void drawBattery() {
    int x = 295, y = 2, w = 20, h = 10;
    uint8_t level = boot_GetBatteryStatus();

    gfx_SetColor(0);
    gfx_Rectangle(x, y, w, h);
    gfx_FillRectangle(x + w, y + 3, 2, 4);

    int fillW = (level * (w - 2)) / 4;
    if (level <= 1) gfx_SetColor(224);
    else if (level <= 2) gfx_SetColor(231);
    else gfx_SetColor(7);

    if (fillW > 0) gfx_FillRectangle(x + 1, y + 1, fillW, h - 2);
}

void drawTimer() {
    clock_t now = clock();
    int elapsed = (now - lastActivityTime) / CLOCKS_PER_SEC;
    int remaining = INACTIVITY_TIMEOUT - elapsed;
    if (remaining < 0) remaining = 0;

    int mins = remaining / 60;
    int secs = remaining % 60;

    char buf[10];
    buf[0] = '0' + mins;
    buf[1] = ':';
    buf[2] = '0' + (secs / 10);
    buf[3] = '0' + (secs % 10);
    buf[4] = '\0';

    gfx_SetTextFGColor(remaining <= 30 ? 224 : 0);
    gfx_PrintStringXY(buf, 5, 3);
}

void drawTable() {
    gfx_FillScreen(255);

    drawBattery();
    drawTimer();

    gfx_SetTextFGColor(0);
    gfx_SetTextScale(1, 1);
    gfx_PrintStringXY("PROJECTILE MOTION", 95, 3);
    
    int startX = 5;
    int startY = 16;
    int colW = 68;
    int rowH = 15;
    int labelW = 20;
    
    gfx_SetColor(0);
    gfx_PrintStringXY("X", startX + labelW + 28, startY + 2);
    gfx_PrintStringXY("Y", startX + labelW + colW + 28, startY + 2);
    
    gfx_HorizLine(startX, startY + rowH - 2, labelW + colW * 2 + 10);
    
    for (int row = 0; row < ROWS; row++) {
        int y = startY + rowH + row * rowH;
        
        gfx_SetTextFGColor(0);
        gfx_PrintStringXY(rowLabels[row], startX + 2, y + 3);
        
        for (int col = 0; col < COLS; col++) {
            int x = startX + labelW + col * colW;
            
            bool selected = (row == curRow && col == curCol && !inputMode);
            bool editing = (row == curRow && col == curCol && inputMode);
            
            if (selected) {
                gfx_SetColor(183);
                gfx_FillRectangle(x, y, colW - 2, rowH - 2);
            } else if (editing) {
                gfx_SetColor(239);
                gfx_FillRectangle(x, y, colW - 2, rowH - 2);
            }
            
            gfx_SetColor(0);
            gfx_Rectangle(x, y, colW - 2, rowH - 2);
            
            char valStr[15];
            float* vals = (col == 0) ? xVals : yVals;
            bool* known = (col == 0) ? xKnown : yKnown;
            bool* userSet = (col == 0) ? xUserSet : yUserSet;
            
            if (editing) {
                gfx_SetTextFGColor(0);
                char displayStr[15];
                if (isNegative && inputLen > 0) {
                    displayStr[0] = '-';
                    strcpy(displayStr + 1, inputBuf);
                } else if (isNegative) {
                    strcpy(displayStr, "-");
                } else if (inputLen == 0) {
                    strcpy(displayStr, "_");
                } else {
                    strcpy(displayStr, inputBuf);
                }
                gfx_PrintStringXY(displayStr, x + 3, y + 3);
            } else if (known[row]) {
                floatToStr(vals[row], valStr);
                gfx_SetTextFGColor(userSet[row] ? 0 : 24);
                gfx_PrintStringXY(valStr, x + 3, y + 3);
            } else {
                gfx_SetTextFGColor(0);
                gfx_PrintStringXY("?", x + 28, y + 3);
            }
        }
    }
    
    int extraStartY = startY + rowH;
    int extraRowH = 15;
    int boxW = 60;
    int rightColX = 168;
    
    const char* extraLabels[3] = {"v0:", "ang:", "vf:"};
    const char* extraUnits[3] = {"m/s", "deg", "m/s"};
    float extraVals[3] = {launchSpeed, launchAngle, finalSpeed};
    bool extraKnown[3] = {speedKnown, angleKnown, finalSpeedKnown};
    bool extraUserSet[3] = {speedUserSet, angleUserSet, finalSpeedUserSet};
    
    for (int i = 0; i < 3; i++) {
        int y = extraStartY + i * extraRowH;
        
        gfx_SetTextFGColor(0);
        gfx_PrintStringXY(extraLabels[i], rightColX, y + 3);
        
        int boxX = rightColX + 28;
        bool selected = (curRow == ROWS + i && !inputMode);
        bool editing = (curRow == ROWS + i && inputMode);
        
        if (selected) {
            gfx_SetColor(183);
            gfx_FillRectangle(boxX, y, boxW, extraRowH - 2);
        } else if (editing) {
            gfx_SetColor(239);
            gfx_FillRectangle(boxX, y, boxW, extraRowH - 2);
        }
        gfx_SetColor(0);
        gfx_Rectangle(boxX, y, boxW, extraRowH - 2);
        
        if (editing) {
            gfx_SetTextFGColor(0);
            char displayStr[15];
            if (isNegative && inputLen > 0) {
                displayStr[0] = '-';
                strcpy(displayStr + 1, inputBuf);
            } else if (isNegative) {
                strcpy(displayStr, "-");
            } else if (inputLen == 0) {
                strcpy(displayStr, "_");
            } else {
                strcpy(displayStr, inputBuf);
            }
            gfx_PrintStringXY(displayStr, boxX + 3, y + 3);
        } else if (extraKnown[i]) {
            char valStr[15];
            floatToStr(extraVals[i], valStr);
            gfx_SetTextFGColor(extraUserSet[i] ? 0 : 24);
            gfx_PrintStringXY(valStr, boxX + 3, y + 3);
        } else {
            gfx_SetTextFGColor(0);
            gfx_PrintStringXY("?", boxX + 25, y + 3);
        }
        
        gfx_SetTextFGColor(0);
        gfx_PrintStringXY(extraUnits[i], boxX + boxW + 3, y + 3);
    }
    
    int bottomY = startY + rowH + ROWS * rowH + 6;
    
    gfx_SetTextFGColor(0);
    
    if (yKnown[2] && yKnown[4]) {
        float v0y = yVals[2];
        float ay = yVals[4];
        float y0 = yKnown[0] ? yVals[0] : 0;
        float maxH;
        if (ay < 0 && v0y > 0) {
            float tMax = -v0y / ay;
            maxH = y0 + v0y * tMax + 0.5f * ay * tMax * tMax;
        } else {
            maxH = y0;
        }
        char buf[20];
        gfx_PrintStringXY("MH:", 168, extraStartY + 3 * extraRowH + 2);
        floatToStr(maxH, buf);
        gfx_PrintStringXY(buf, 192, extraStartY + 3 * extraRowH + 2);
        gfx_PrintStringXY("m", 260, extraStartY + 3 * extraRowH + 2);
    }
    if (xKnown[6]) {
        char buf[20];
        gfx_PrintStringXY("ToF:", 168, extraStartY + 4 * extraRowH + 2);
        floatToStr(xVals[6], buf);
        gfx_PrintStringXY(buf, 200, extraStartY + 4 * extraRowH + 2);
        gfx_PrintStringXY("s", 260, extraStartY + 4 * extraRowH + 2);
    }
    
    gfx_SetTextFGColor(24);
    char allEqs[256] = "";
    if (xEqUsed[0]) strcpy(allEqs, xEqUsed);
    if (yEqUsed[0]) {
        char* tok = yEqUsed;
        while (*tok) {
            char* end = strchr(tok, '|');
            char eq[64];
            if (end) {
                int len = end - tok;
                strncpy(eq, tok, len);
                eq[len] = '\0';
                tok = end + 1;
            } else {
                strcpy(eq, tok);
                tok += strlen(tok);
            }
            if (strstr(allEqs, eq) == NULL) {
                if (allEqs[0]) strcat(allEqs, "|");
                strcat(allEqs, eq);
            }
        }
    }

    int eqY = bottomY;
    char* tok = allEqs;
    while (*tok) {
        char* end = strchr(tok, '|');
        char eq[64];
        if (end) {
            int len = end - tok;
            strncpy(eq, tok, len);
            eq[len] = '\0';
            tok = end + 1;
        } else {
            strcpy(eq, tok);
            tok += strlen(tok);
        }
        gfx_PrintStringXY(eq, 5, eqY);
        eqY += 10;
    }
    
    gfx_SetTextFGColor(160);
    gfx_PrintStringXY("Built: " __DATE__ " " __TIME__, 5, 230);
    
    int miniX = 165;
    int miniY = 108;
    int miniW = 150;
    int miniH = 60;
    
    gfx_SetColor(200);
    gfx_Rectangle(miniX, miniY, miniW, miniH);
    
    if (xKnown[2] && yKnown[2] && xKnown[6] && xVals[6] > 0) {
        float totalTime = xVals[6];
        float x0 = xKnown[0] ? xVals[0] : 0;
        float y0 = yKnown[0] ? yVals[0] : 0;
        float maxPx = x0, minPx = x0, maxPy = y0, minPy = y0;
        for (int i = 0; i <= 20; i++) {
            float tt = (i / 20.0f) * totalTime;
            float px = x0 + xVals[2] * tt + 0.5f * xVals[4] * tt * tt;
            float py = y0 + yVals[2] * tt + 0.5f * yVals[4] * tt * tt;
            if (px > maxPx) maxPx = px;
            if (px < minPx) minPx = px;
            if (py > maxPy) maxPy = py;
            if (py < minPy) minPy = py;
        }
        float rangeX = maxPx - minPx;
        float rangeY = maxPy - minPy;
        if (rangeX < 1) rangeX = 1;
        if (rangeY < 1) rangeY = 1;
        minPx -= rangeX * 0.1f; maxPx += rangeX * 0.1f;
        minPy -= rangeY * 0.1f; maxPy += rangeY * 0.1f;
        rangeX = maxPx - minPx;
        rangeY = maxPy - minPy;
        
        gfx_SetColor(24);
        int lastSx = -1, lastSy = -1;
        for (int i = 0; i <= 30; i++) {
            float tt = (i / 30.0f) * totalTime;
            float px = x0 + xVals[2] * tt + 0.5f * xVals[4] * tt * tt;
            float py = y0 + yVals[2] * tt + 0.5f * yVals[4] * tt * tt;
            int sx = miniX + 5 + (int)(((px - minPx) / rangeX) * (miniW - 10));
            int sy = miniY + miniH - 5 - (int)(((py - minPy) / rangeY) * (miniH - 10));
            if (lastSx >= 0) gfx_Line(lastSx, lastSy, sx, sy);
            lastSx = sx; lastSy = sy;
        }
        gfx_SetColor(224);
        int startSx = miniX + 5 + (int)(((x0 - minPx) / rangeX) * (miniW - 10));
        int startSy = miniY + miniH - 5 - (int)(((y0 - minPy) / rangeY) * (miniH - 10));
        gfx_FillCircle(startSx, startSy, 2);
    }
    
    gfx_SetTextFGColor(0);
    gfx_PrintStringXY("p0\\pf: init\\final pos (m)", 140, 172);
    gfx_PrintStringXY("v0\\vf: init\\final vel (m/s)", 140, 181);
    gfx_PrintStringXY("a: acceleration (m/s^2)", 140, 190);
    gfx_PrintStringXY("d: displacement (m)", 140, 199);
    gfx_PrintStringXY("MH: maximum height (m)", 140, 208);
    gfx_PrintStringXY("ToF: time of flight (s)", 140, 217);
    gfx_PrintStringXY("[graph]", 265, 230);
}

void drawGraph() {
    if (!xKnown[2] || !yKnown[2] || !xKnown[6]) return;
    
    gfx_FillScreen(255);
    gfx_SetColor(0);
    
    float totalTime = xVals[6];
    if (totalTime <= 0) totalTime = 5.0f;
    
    float x0 = xKnown[0] ? xVals[0] : 0;
    float y0 = yKnown[0] ? yVals[0] : 0;
    float maxX = x0, minX = x0, maxY = y0, minY = y0;
    
    for (int i = 0; i <= 50; i++) {
        float t = (i / 50.0f) * totalTime;
        float px = x0 + xVals[2] * t + 0.5f * xVals[4] * t * t;
        float py = y0 + yVals[2] * t + 0.5f * yVals[4] * t * t;
        if (px > maxX) maxX = px;
        if (px < minX) minX = px;
        if (py > maxY) maxY = py;
        if (py < minY) minY = py;
    }
    
    float rangeX = maxX - minX;
    float rangeY = maxY - minY;
    if (rangeX < 1) rangeX = 1;
    if (rangeY < 1) rangeY = 1;
    
    float margin = 0.1f;
    minX -= rangeX * margin;
    maxX += rangeX * margin;
    minY -= rangeY * margin;
    maxY += rangeY * margin;
    rangeX = maxX - minX;
    rangeY = maxY - minY;
    
    int graphX = 30;
    int graphY = 15;
    int graphW = 280;
    int graphH = 190;
    
    gfx_SetColor(200);
    gfx_Rectangle(graphX, graphY, graphW, graphH);
    
    gfx_SetColor(0);
    int zeroX = graphX + (int)((-minX / rangeX) * graphW);
    int zeroY = graphY + graphH - (int)((-minY / rangeY) * graphH);
    
    if (zeroX >= graphX && zeroX <= graphX + graphW) {
        gfx_VertLine(zeroX, graphY, graphH);
    }
    if (zeroY >= graphY && zeroY <= graphY + graphH) {
        gfx_HorizLine(graphX, zeroY, graphW);
    }
    
    gfx_SetTextFGColor(0);
    gfx_PrintStringXY("x(m)", 290, graphY + graphH + 5);
    gfx_PrintStringXY("y", graphX - 15, graphY);
    
    gfx_SetColor(24);
    
    int lastSx = -1, lastSy = -1;
    for (int i = 0; i <= 100; i++) {
        float t = (i / 100.0f) * totalTime;
        float px = x0 + xVals[2] * t + 0.5f * xVals[4] * t * t;
        float py = y0 + yVals[2] * t + 0.5f * yVals[4] * t * t;
        
        int sx = graphX + (int)(((px - minX) / rangeX) * graphW);
        int sy = graphY + graphH - (int)(((py - minY) / rangeY) * graphH);
        
        if (sx >= graphX && sx <= graphX + graphW && sy >= graphY && sy <= graphY + graphH) {
            if (lastSx >= 0) {
                gfx_Line(lastSx, lastSy, sx, sy);
            }
            lastSx = sx;
            lastSy = sy;
        }
    }
    
    gfx_SetColor(224);
    int startSx = graphX + (int)(((x0 - minX) / rangeX) * graphW);
    int startSy = graphY + graphH - (int)(((y0 - minY) / rangeY) * graphH);
    gfx_FillCircle(startSx, startSy, 4);
    
    gfx_SetTextFGColor(0);
    gfx_PrintStringXY("Any key to return", 110, 225);
    
    gfx_BlitBuffer();
    
    while (!kb_AnyKey()) kb_Scan();
    while (kb_AnyKey()) kb_Scan();
}

void startInput() {
    inputMode = true;
    inputLen = 0;
    inputBuf[0] = '\0';
    hasDecimal = false;
    isNegative = false;
}

void finishInput() {
    if (inputLen > 0 || isNegative) {
        char fullStr[15];
        if (isNegative) {
            fullStr[0] = '-';
            strcpy(fullStr + 1, inputBuf);
        } else {
            strcpy(fullStr, inputBuf);
        }
        
        char* end;
        real_t r = os_StrToReal(fullStr, &end);
        float val = os_RealToFloat(&r);
        
        if (curRow < ROWS) {
            if (curCol == 0) {
                xVals[curRow] = val;
                xKnown[curRow] = true;
                xUserSet[curRow] = true;
            } else {
                yVals[curRow] = val;
                yKnown[curRow] = true;
                yUserSet[curRow] = true;
            }
        } else if (curRow == ROWS) {
            launchSpeed = val;
            speedKnown = true;
            speedUserSet = true;
        } else if (curRow == ROWS + 1) {
            launchAngle = val;
            angleKnown = true;
            angleUserSet = true;
        } else {
            finalSpeed = val;
            finalSpeedKnown = true;
            finalSpeedUserSet = true;
        }
        
        autoSolve();
    }
    inputMode = false;
}

void cancelInput() {
    inputMode = false;
}

void clearCell() {
    if (curRow < ROWS) {
        if (curCol == 0) {
            xKnown[curRow] = false;
            xUserSet[curRow] = false;
            xVals[curRow] = 0;
        } else {
            yKnown[curRow] = false;
            yUserSet[curRow] = false;
            yVals[curRow] = 0;
        }
    } else if (curRow == ROWS) {
        speedKnown = false;
        speedUserSet = false;
        launchSpeed = 0;
    } else if (curRow == ROWS + 1) {
        angleKnown = false;
        angleUserSet = false;
        launchAngle = 0;
    } else {
        finalSpeedKnown = false;
        finalSpeedUserSet = false;
        finalSpeed = 0;
    }
    autoSolve();
}

void resetAll() {
    initData();
}

int main(void) {
    gfx_Begin();
    gfx_SetDrawBuffer();

    initData();
    lastActivityTime = clock();

    bool running = true;
    bool prevUp = false, prevDown = false, prevLeft = false, prevRight = false;
    bool prevEnter = false, prevClear = false, prevDel = false;
    bool prevMode = false, prevGraph = false;
    bool prevKeys[10] = {false};
    bool prevNeg = false, prevDot = false;

    while (running) {
        drawTable();
        gfx_BlitBuffer();

        kb_Scan();

        bool up = kb_Data[7] & kb_Up;
        bool down = kb_Data[7] & kb_Down;
        bool left = kb_Data[7] & kb_Left;
        bool right = kb_Data[7] & kb_Right;
        bool enter = kb_Data[6] & kb_Enter;
        bool clear = kb_Data[6] & kb_Clear;
        bool del = kb_Data[1] & kb_Del;
        bool mode = kb_Data[1] & kb_Mode;
        bool graph = kb_Data[1] & kb_Graph;

        bool keys[10];
        keys[0] = kb_Data[3] & kb_0;
        keys[1] = kb_Data[3] & kb_1;
        keys[2] = kb_Data[4] & kb_2;
        keys[3] = kb_Data[5] & kb_3;
        keys[4] = kb_Data[3] & kb_4;
        keys[5] = kb_Data[4] & kb_5;
        keys[6] = kb_Data[5] & kb_6;
        keys[7] = kb_Data[3] & kb_7;
        keys[8] = kb_Data[4] & kb_8;
        keys[9] = kb_Data[5] & kb_9;
        bool neg = kb_Data[5] & kb_Chs;
        bool dot = kb_Data[4] & kb_DecPnt;

        if (kb_AnyKey()) lastActivityTime = clock();

        clock_t now = clock();
        if ((now - lastActivityTime) / CLOCKS_PER_SEC >= INACTIVITY_TIMEOUT) {
            running = false;
            continue;
        }
        
        if (!inputMode) {
            if (up && !prevUp) {
                if (curRow >= ROWS) {
                    if (curRow == ROWS) curRow = ROWS + 2;
                    else curRow--;
                } else {
                    curRow = (curRow - 1 + ROWS) % ROWS;
                }
            }
            if (down && !prevDown) {
                if (curRow >= ROWS) {
                    if (curRow == ROWS + 2) curRow = ROWS;
                    else curRow++;
                } else {
                    curRow = (curRow + 1) % ROWS;
                }
            }
            if (left && !prevLeft) {
                if (curRow >= ROWS) {
                    curRow = curRow - ROWS;
                    curCol = 1;
                } else if (curRow < 3 && curCol == 0) {
                    curRow = ROWS + curRow;
                } else {
                    curCol = (curCol - 1 + COLS) % COLS;
                }
            }
            if (right && !prevRight) {
                if (curRow < 3 && curCol == 1) {
                    curRow = ROWS + curRow;
                } else if (curRow >= ROWS) {
                    curRow = curRow - ROWS;
                    curCol = 0;
                } else {
                    curCol = (curCol + 1) % COLS;
                }
            }
            
            if (enter && !prevEnter) startInput();
            if (del && !prevDel) clearCell();
            if (mode && !prevMode) resetAll();
            if (graph && !prevGraph) drawGraph();
            if (clear && !prevClear) running = false;
            
            for (int i = 0; i <= 9; i++) {
                if (keys[i] && !prevKeys[i]) {
                    startInput();
                    inputBuf[0] = '0' + i;
                    inputBuf[1] = '\0';
                    inputLen = 1;
                    break;
                }
            }
            if (neg && !prevNeg) {
                startInput();
                isNegative = true;
            }
            if (dot && !prevDot) {
                startInput();
                inputBuf[0] = '.';
                inputBuf[1] = '\0';
                inputLen = 1;
                hasDecimal = true;
            }
        } else {
            for (int i = 0; i <= 9; i++) {
                if (keys[i] && !prevKeys[i] && inputLen < 10) {
                    inputBuf[inputLen++] = '0' + i;
                    inputBuf[inputLen] = '\0';
                }
            }
            
            if (dot && !prevDot && !hasDecimal && inputLen < 10) {
                inputBuf[inputLen++] = '.';
                inputBuf[inputLen] = '\0';
                hasDecimal = true;
            }
            
            if (neg && !prevNeg) {
                isNegative = !isNegative;
            }
            
            if (del && !prevDel && inputLen > 0) {
                inputLen--;
                if (inputBuf[inputLen] == '.') hasDecimal = false;
                inputBuf[inputLen] = '\0';
            }
            
            if (enter && !prevEnter) finishInput();
            if (clear && !prevClear) cancelInput();
        }
        
        prevUp = up; prevDown = down; prevLeft = left; prevRight = right;
        prevEnter = enter; prevClear = clear; prevDel = del;
        prevMode = mode; prevGraph = graph;
        for (int i = 0; i <= 9; i++) prevKeys[i] = keys[i];
        prevNeg = neg; prevDot = dot;
    }
    
    gfx_End();
    return 0;
}
