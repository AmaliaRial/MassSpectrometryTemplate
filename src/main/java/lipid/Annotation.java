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
    private final IonizationMode ionizationMode;
    private String adduct;
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode,Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
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

    public IonizationMode getIonizationMode() { return  ionizationMode; }

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

    /**
     * Automatically detects the most probable adduct comparing the difference in mass
     * between the grouped peaks (groupedSignals) with the characteristic masses of the known adducts.
     *
     * @param ppmTolerance tolerance in parts per million to be able to consider a valid match
     */
    public void detectAdduct(int ppmTolerance) {
        if (groupedSignals.size() < 2) {
            this.adduct = null;
            return;
        }

        List<Peak> peaks = new ArrayList<>(groupedSignals);
        peaks.sort(Comparator.comparingDouble(Peak::getMz)); // order by m/z

        // combines both maps of adducts
        Map<String, Double> adducts = new LinkedHashMap<>();
        if (ionizationMode == IonizationMode.POSITIVE) {
            adducts.putAll(AdductList.MAPMZPOSITIVEADDUCTS);
        } else {
            adducts.putAll(AdductList.MAPMZNEGATIVEADDUCTS);
        }

        // By default we assume that the adduct with the lowest m/z is the main one
        Peak basePeak = peaks.get(0);

        String bestMatch = null;

        for (int i = 0; i < peaks.size(); i++) {
            for (int j = i + 1; j < peaks.size(); j++) {
                Peak p1 = peaks.get(i);
                Peak p2 = peaks.get(j);
                double deltaMz = Math.abs(p2.getMz() - p1.getMz());

                // We compare deltaMz against the difference between all pairs of known adducts
                for (String a1 : adducts.keySet()) {
                    for (String a2 : adducts.keySet()) {
                        if (a1.equals(a2)) continue;

                        double mass1 = Math.abs(adducts.get(a1));
                        double mass2 = Math.abs(adducts.get(a2));
                        double expectedDiff = Math.abs(mass1 - mass2);

                        int ppm = Adduct.calculatePPMIncrement(deltaMz, expectedDiff);
                        if (ppm <= ppmTolerance) {
                            // We assume that the one with the lowest m/z is the one we will use as a reference
                            bestMatch = p1.getMz() <= p2.getMz() ? a1 : a2;
                            this.adduct = bestMatch;
                            return;
                        }
                    }
                }
            }
        }

        // Si no se detectó nada
        this.adduct = null;
    }


}
