package splitset

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"regexp"

	"github.com/Open-Science-Global/poly/transform"
	"github.com/schwarmco/go-cartesian-product"
)

type Overhang struct {
	Sequence string
	Start    int
	End      int
}

type Region struct {
	Sequence string
	Start    int
	End      int
}

type Permutation struct {
	Overhangs  []string
	Matches    int
	Mismatches int
}

type Fragment struct {
	StartOverhang    string
	EndOverhang      string
	SequenceFragment string
}

func SplitSet(sequence string, numberOfFragments int, spaceAround int) []Fragment {
	regions := GetRegions(sequence, numberOfFragments, spaceAround)
	overhangs := FindBestOverhangs(regions)
	positions := FindPositions(regions, overhangs)

	return GetFragments(sequence, overhangs, positions)
}

func GetFragments(sequence string, permutation Permutation, positions []Overhang) []Fragment {
	startPoint := 0
	lastOverhang := ""
	var fragments []Fragment

	for _, overhang := range permutation.Overhangs {

		for _, position := range positions {
			if overhang == position.Sequence {
				if lastOverhang == "" {
					lastOverhang = overhang
					fragments = append(fragments, Fragment{"", lastOverhang, sequence[startPoint:position.End]})
					startPoint = position.Start
				} else {
					fragments = append(fragments, Fragment{lastOverhang, overhang, sequence[startPoint:position.End]})
					startPoint = position.Start
					lastOverhang = overhang
				}
				break
			}
		}
	}

	fragments = append(fragments, Fragment{lastOverhang, "", sequence[startPoint:]})

	return fragments
}

func FindPositions(regions []Region, overhangs Permutation) []Overhang {
	var locations []Overhang
	for index, region := range regions {
		overhang := overhangs.Overhangs[index]
		re := regexp.MustCompile(overhang)
		locs := re.FindAllStringIndex(region.Sequence, -1)
		for _, loc := range locs {
			start := loc[0]
			end := loc[1]
			locations = append(locations, Overhang{region.Sequence[start:end], region.Start + start, region.Start + end})
		}
	}
	return locations
}

func GetRegions(sequence string, numberOfFragments int, spaceAround int) []Region {
	fragmentSize := len(sequence) / numberOfFragments

	var regions []Region
	for i := 1; i < numberOfFragments; i++ {
		middle := fragmentSize * i
		start := middle - spaceAround
		end := middle + spaceAround
		regions = append(regions, Region{sequence[start:end], start, end})
	}
	return regions
}

func FindBestOverhangs(regions []Region) Permutation {
	freqTable := ReadFreqTableJSON("./dataset/freq_overhang.json")
	permutations := generatePermutations(regions)
	var infos []Permutation
	for _, permutation := range permutations {
		info := getInfo(permutation, freqTable)
		infos = append(infos, info)
	}

	bestPermutation := filterBestPermutationByInfo(infos)

	return bestPermutation
}

func filterBestPermutationByInfo(infos []Permutation) Permutation {
	bestPermutation := Permutation{}
	matches := 0
	mismatches := 100
	for _, info := range infos {
		if matches <= info.Matches && mismatches >= info.Mismatches {
			bestPermutation = info
			matches = info.Matches
			mismatches = info.Mismatches
			fmt.Println(bestPermutation)
		}
	}
	return bestPermutation
}

func generatePermutations(regions []Region) [][]string {
	var kmerTables [][]interface{}
	for _, region := range regions {
		kmerTables = append(kmerTables, kmerToList(transform.GetKmerTable(4, region.Sequence)))
	}

	products := cartesian.Iter(kmerTables...)
	var permutations [][]string
	for product := range products {
		var permutation []string
		for _, element := range product {
			permutation = append(permutation, element.(string))
		}
		if uniquePermutation(permutation) {
			permutations = append(permutations, permutation)
		}
	}

	return permutations
}

func uniquePermutation(permutation []string) bool {
	set := make(map[string]bool)
	for _, overhang := range permutation {
		_, exist := set[overhang]
		if exist {
			return false
		}
		set[overhang] = true
		set[transform.ReverseComplement(overhang)] = true
	}
	return true
}

func kmerToList(kmerTable map[string]bool) []interface{} {
	var list []interface{}
	for key, _ := range kmerTable {
		list = append(list, key)
	}
	return list
}

func getInfo(options []string, freqTable map[string]map[string]int) Permutation {
	var optionsWithReverseComplement []string
	matches := 0
	mismatches := 0

	for _, option := range options {
		optionsWithReverseComplement = append(optionsWithReverseComplement, option, transform.ReverseComplement(option))
	}

	for _, topStrand := range optionsWithReverseComplement {
		for _, bottomStrand := range optionsWithReverseComplement {
			if transform.ReverseComplement(topStrand) == bottomStrand && topStrand != transform.ReverseComplement(topStrand) {
				matches += freqTable[topStrand][bottomStrand]
			} else {
				mismatches += freqTable[topStrand][bottomStrand]
			}
		}
	}
	return Permutation{options, matches, mismatches}
}

func ParseFreqTableJSON(file []byte) map[string]map[string]int {
	var codontable map[string]map[string]int
	_ = json.Unmarshal([]byte(file), &codontable)
	return codontable
}

// ReadCodonJSON reads a Table JSON file.
func ReadFreqTableJSON(path string) map[string]map[string]int {
	file, _ := ioutil.ReadFile(path)
	codontable := ParseFreqTableJSON(file)
	return codontable
}
