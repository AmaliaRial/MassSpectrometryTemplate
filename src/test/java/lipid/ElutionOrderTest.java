package lipid;
import org.drools.ruleunits.api.RuleUnitInstance;
import org.drools.ruleunits.api.RuleUnitProvider;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class ElutionOrderTest {


    static final Logger LOG = LoggerFactory.getLogger(ElutionOrderTest.class);

    // !!TODO For the adduct detection both regular algorithms or drools can be used as far the tests are passed.


    @Before
    public void setup() {
        // !! TODO Empty by now,you can create common objects for all tests.
    }


    @Test
    public void score1BasedOnRT() {
        // Assume lipids already annotated
        LOG.info("Creating RuleUnit");
        LipidScoreUnit lipidScoreUnit = new LipidScoreUnit();

        RuleUnitInstance<LipidScoreUnit> instance = RuleUnitProvider.get().createRuleUnitInstance(lipidScoreUnit);


        // TODO CHECK THE Monoisotopic MASSES OF THE COMPOUNDS IN https://chemcalc.org/
        Lipid lipid1 = new Lipid(1, "TG 54:3", "C57H104O6", "TG", 54, 3); // MZ of [M+H]+ = 885.79057
        Lipid lipid2 = new Lipid(2, "TG 52:3", "C55H100O6", "TG", 52, 3); // MZ of [M+H]+ = 857.75927
        Lipid lipid3 = new Lipid(3, "TG 56:3", "C59H108O6", "TG", 56, 3); // MZ of [M+H]+ = 913.82187
        Annotation annotation1 = new Annotation(lipid1, 885.79056, 10E6, 10d);
        Annotation annotation2 = new Annotation(lipid2, 857.7593, 10E7, 9d);
        Annotation annotation3 = new Annotation(lipid3, 913.822, 10E5, 11d);


        LOG.info("Insert data");

        try {
            lipidScoreUnit.getAnnotations().add(annotation1);
            lipidScoreUnit.getAnnotations().add(annotation2);
            lipidScoreUnit.getAnnotations().add(annotation3);

            LOG.info("Run query. Rules are also fired");
            instance.fire();

            // Here the logic that we expect. In this case we expect the full 3 annotations to have a positive score of 1

            assertEquals(1.0, annotation1.getNormalizedScore(), 0.01);
            assertEquals(1.0, annotation2.getNormalizedScore(), 0.01);
            assertEquals(1.0, annotation3.getNormalizedScore(), 0.01);

        }
        finally {
            instance.close();
        }
    }

}
