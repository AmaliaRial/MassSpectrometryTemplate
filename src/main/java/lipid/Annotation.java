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

        // If there are less than two signals, no difference can be calculated → no adduct can be inferred
        if (groupedSignals.size() < 2) {
            this.adduct = null;
            return;
        }

        // We convert the set of signals into a list ordered by mass/charge (mz)
        List<Peak> peaks = new ArrayList<>(groupedSignals);
        peaks.sort(Comparator.comparingDouble(Peak::getMz)); // sorts the peaks by their m/z value

        // We compare all possible combinations of two different peaks
        for (int i = 0; i < peaks.size(); i++) {
            for (int j = i + 1; j < peaks.size(); j++) {
                Peak lower = peaks.get(i);           // peak with lower mz
                Peak higher = peaks.get(j);          // peak with higher mz
                double deltaMz = higher.getMz() - lower.getMz(); // difference in m/z

                // Run through the maps of known positive and negative adducts
                for (Map<String, Double> map : List.of(
                        AdductList.MAPMZPOSITIVEADDUCTS,
                        AdductList.MAPMZNEGATIVEADDUCTS)) {

                    // For each adduct, we get its name and associated mass
                    for (Map.Entry<String, Double> entry : map.entrySet()) {
                        String candidate = entry.getKey();                  // adduct name
                        double adductMassDiff = Math.abs(entry.getValue()); // difference in mass

                        // We calculate the difference between deltaMz and the expected mass, in ppm
                        int ppm = Adduct.calculatePPMIncrement(deltaMz, adductMassDiff);

                        // If it is within the allowed tolerance, we consider it a match
                        if (ppm <= ppmTolerance) {
                            this.adduct = candidate; // assign the detected adduct
                            return;
                        }
                    }
                }
            }
        }

        // If no combination matches a known adduct, we set it to null
        this.adduct = null;
    }


}
