package splitset

import "fmt"

func ExampleSplitSet() {
	result := SplitSet("TGGTACGAAAATTAGGGGATCTACCTAGAAAGCCACAAGGCGATAGGTCAAGCTTAAAGAACCCTTACATGGATCTTACAGATTCTGAAAGTAAAGAAACAACAGAGGTTAAACAAACAGAACCAAAAAGAAAAAAAGCATTGTTGAAAACAATGAAAGTTGATGTTTCAATCCATAATAAGATTAAATCGCTGCACGAAATTCTGGCAGCATCCGAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3, 25)
	fmt.Println(result)
	//Output: [CATG TGTT]
}

func ExampleReadFreqTableJSON() {
	result := ReadFreqTableJSON("./dataset/freq_overhang_clean.json")
	fmt.Println(result["TCTA"]["TAGA"])

	// Output: 210
}

func ExampleGeneratePermutations() {
	regions := []Region{Region{"AAAAT", 0, 0}, Region{"TTTT", 0, 0}, Region{"CCCC", 0, 0}}
	result := generatePermutations(regions)

	fmt.Println(result)
	// Output: [[AAAA TTTT CCCC] [AAAT TTTT CCCC]]
}

// Check with Isaac Larkin, this test case is actually pretty good
func ExampleFindBestOverhangs() {
	regions := []Region{Region{"AAAAT", 0, 0}, Region{"TTTT", 0, 0}, Region{"CCCC", 0, 0}}
	permutations := findBestOverhangs(regions)

	fmt.Println(permutations)
	// Output: [AAAA TTTT CCCC]
}

func ExampleFindBestOverhangsReal() {
	regions := getRegions("TGGTACGAAAATTAGGGGATCTACCTAGAAAGCCACAAGGCGATAGGTCAAGCTTAAAGAACCCTTACATGGATCTTACAGATTCTGAAAGTAAAGAAACAACAGAGGTTAAACAAACAGAACCAAAAAGAAAAAAAGCATTGTTGAAAACAATGAAAGTTGATGTTTCAATCCATAATAAGATTAAATCGCTGCACGAAATTCTGGCAGCATCCGAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3, 25)
	permutations := findBestOverhangs(regions)

	fmt.Println(permutations)
	// Output: [AAAA TTTT CCCC]
}

func ExampleGetRegion() {
	result := getRegions("TGGTACGAAAATTAGGGGATCTACCTAGAAAGCCACAAGGCGATAGGTCAAGCTTAAAGAACCCTTACATGGATCTTACAGATTCTGAAAGTAAAGAAACAACAGAGGTTAAACAAACAGAACCAAAAAGAAAAAAAGCATTGTTGAAAACAATGAAAGTTGATGTTTCAATCCATAATAAGATTAAATCGCTGCACGAAATTCTGGCAGCATCCGAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAA", 3, 25)

	fmt.Println(result)
	// Output: [AGAACCCTTACATGGATCTTACAGATTCTGAAAGTAAAGAAACAACAGAG CATTGTTGAAAACAATGAAAGTTGATGTTTCAATCCATAATAAGATTAAA]
}

func ExampleUniquePermutation() {
	result := uniquePermutation([]string{"AACA", "TGTT"})
	fmt.Println(result)
	// Output: false
}

func ExampleUniquePermutationDuplicate() {
	result := uniquePermutation([]string{"AAAA", "TTTT", "CCCC"})
	fmt.Println(result)
	// Output: false
}
