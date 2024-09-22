
using UdonSharp;
using UnityEngine;
using VRC.SDKBase;
using VRC.Udon;
using System;

using static UnityEngine.Mathf;

public class moon_movement : UdonSharpBehaviour
{
    public float DistanceMultiplier = 1;

    Vector3 moon_position(float days_since_J2000) {

        // ECI position (meters) of the Moon at the given MJD (UTC).
        // ported from public domain code:
        // https://possiblywrong.wordpress.com/2014/01/04/computing-positions-and-eclipses-of-the-sun-and-moon/
        // Reference: Montenbruck and Gill, Satellite Orbits. Berlin: Springer,
        // 2005, Chapter 3.3.2.
        // float days_since_J2000 = mjd_utc - 51544.5 + dt_tt / 86400;
        float T = days_since_J2000 / 36525.0f; // centuries since J2000

        // equation (3.47)
        float L0 = Deg2Rad*(218.31617f + 481267.88088f * T);
        float l = Deg2Rad*(134.96292f + 477198.86753f * T);
        float lp = Deg2Rad*(357.52543f + 35999.04944f * T);
        float F = Deg2Rad*(93.27283f + 483202.01873f * T);
        float D = Deg2Rad*(297.85027f + 445267.11135f * T);

        // equation (3.48)
        float dL = Deg2Rad*((22640f * Sin(l) + 769f * Sin(2 * l)
            - 4586f * Sin(l - 2f * D) + 2370f * Sin(2f * D)
            - 668f * Sin(lp) - 412 * Sin(2f * F)
            - 212f * Sin(2 * l - 2 * D) - 206 * Sin(l + lp - 2 * D)
            + 192f * Sin(l + 2 * D) - 165 * Sin(lp - 2 * D)
            + 148f * Sin(l - lp) - 125 * Sin(D)
            - 110f * Sin(l + lp) - 55 * Sin(2 * F - 2 * D)) / 3600);
        float lon = L0 + dL;

        // equation (3.49)
        float beta = Deg2Rad*((18520 * Sin(F + dL + Deg2Rad*(
                (412 * Sin(2 * F) + 541 * Sin(lp)) / 3600))
            - 526 * Sin(F - 2 * D) + 44 * Sin(l + F - 2 * D)
            - 31 * Sin(-l + F - 2 * D) - 25 * Sin(-2 * l + F)
            - 23 * Sin(lp + F - 2 * D) + 21 * Sin(-l + F)
            + 11 * Sin(-lp + F - 2 * D)) / 3600);

        // equation (3.50)
        float r = 1000 * (385000 - 20905 * Cos(l) - 3699 * Cos(2 * D - l)
             - 2956 * Cos(2 * D) - 570 * Cos(2 * l) + 246 * Cos(2 * l - 2 * D)
             - 205 * Cos(lp - 2 * D) - 171 * Cos(l + 2 * D)
             - 152 * Cos(l + lp - 2 * D));

        // equation (3.51)
        float x = r * Cos(lon) * Cos(beta);
        float y = r * Sin(lon) * Cos(beta);
        float z = r * Sin(beta);

        // equation (3.45)
        float epsilon = Deg2Rad*(23.43929111f);
        float s = Sin(-epsilon);
        float c = Cos(-epsilon);
        return new Vector3(
            x,
             y * c + z * s,
            -y * s + z * c
        );
    }
    void Update()
    {
        // Jan 1, 1970 = 621355968000000000.0 ticks.
        double utcSecondsUnix = DateTime.UtcNow.Ticks / 10000000.0 - 62135596800.0;
        double J2000_in_unix_seconds = 946684800.0;
        float days_since_J2000 = (float)((utcSecondsUnix - J2000_in_unix_seconds) / 86400.0);
        Vector3 pos_meters = moon_position(days_since_J2000);
        gameObject.transform.localPosition = DistanceMultiplier * pos_meters / 6371000.0f;
    }
}
