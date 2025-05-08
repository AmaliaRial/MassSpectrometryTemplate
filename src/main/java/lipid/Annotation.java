package lipid;

import adduct.Adduct;
import adduct.AdductList;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private String adduct;
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime) {
        this(lipid, mz, intensity, retentionTime, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.groupedSignals = new TreeSet<>(Comparator.comparingDouble(Peak::getMz));
        this.groupedSignals.addAll(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    public double getNormalizedScore() {
        if (totalScoresApplied == 0) {
            return 0.0;
        }
        double raw = (double) score / totalScoresApplied;
        return Math.min(1.0, Math.max(0.0,raw));}


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    // TENEMOS UN MONTON DE SEÑALES Y VEMOS A VER A CUAL DE LOS ADUCTOS DE LOS CUALES SABEMOS SU MASA SE APROXIMA MAS
    // PARA ASOCIAR ESA SEÑAL A ESE ADUCTO.
    public void detectAdduct(int ppmTolerance) {
        if (groupedSignals.isEmpty()) {
            this.adduct = null;
            return;
        }
        Peak basePeak = groupedSignals.iterator().next();
        for (Map<String, Double> map : List.of(
                AdductList.MAPMZPOSITIVEADDUCTS,
                AdductList.MAPMZNEGATIVEADDUCTS)) {

            for (String candidate : map.keySet()) {
                // calculamos masa monoisotópica según hipótesis
                double mono = Adduct.getMonoisotopicMassFromMZ(basePeak.getMz(), candidate);
                // ppm entre experimental (this.mz) y teórico
                int ppm = Adduct.calculatePPMIncrement(this.mz, mono);
                if (ppm <= ppmTolerance) {
                    this.adduct = candidate;
                    return;
                }
            }
        }
        // si no se encuentra nada:
        this.adduct=null;
    }

}
