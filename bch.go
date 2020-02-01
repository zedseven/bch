// Package bch facilitates Bose-Chaudhuri-Hocquenghem (BCH) codes and error checking.
package bch

import (
	"fmt"
	"math"
)

// EncodingConfig stores all the computed values from encoding that are necessary for decoding later.
type EncodingConfig struct {
	MaxCorrectableErrors int
	N                    int
	Field                GField
}

// GField stores a Galois field for passing between the Encode() and Decode() operations.
type GField struct {
	AlphaTo [1048576]int
	IndexOf [1048576]int
}


func readP(m int) (n int, p [21]uint8, err error) {
	/*fmt.Println("First, enter a value of m such that the code length is")
	fmt.Println("2**(m-1) - 1 < length <= 2**m - 1")
	for do := false; !do || m <= 1 || m >= 21; do = true {
		fmt.Println("Input m (between 2 and 20): ")
		fmt.Scanf("%d", &m)
	}*/
	if m < 2 || m > 20 {
		// TODO: Throw a proper error
		return
	}

	for i := 1; i < m; i++ {
		p[i] = 0
	}
	p[0], p[m] = 1, 1

	switch m {
	case 2, 3, 4, 6, 7, 15:
		p[1] = 1
	case 5, 11:
		p[2] = 1
	case 8:
		p[4], p[5], p[6] = 1, 1, 1
	case 9:
		p[4] = 1
	case 10, 17, 20:
		p[3] = 1
	case 12:
		p[3], p[4], p[7] = 1, 1, 1
	case 13:
		p[1], p[3], p[4] = 1, 1, 1
	case 14:
		p[1], p[11], p[12] = 1, 1, 1
	case 16:
		p[2], p[3], p[5] = 1, 1, 1
	case 18:
		p[7] = 1
	case 19:
		p[1], p[5], p[6] = 1, 1, 1
	}
	fmt.Printf("p(x) = ")
	n = 1
	for i := 0; i <= m; i++ {
		n *= 2
		fmt.Printf("%1d", p[i])
	}
	fmt.Printf("\n")
	n = n / 2 - 1
	/*ninf := (n + 1) / 2 - 1
	for do := false; !do || length > n || length <= ninf; do = true {
		fmt.Printf("Enter the code length (%d < length <= %d): ", ninf, n)
		fmt.Scanf("%d", &length)
	}*/

	return
}

func generateGF(m, n int, p [21]uint8) (alphaTo, indexOf [1048576]int, err error) {
	mask := 1
	alphaTo[m] = 0
	for i := 0; i < m; i++ {
		alphaTo[i] = mask
		indexOf[alphaTo[i]] = i
		if p[i] != 0 {
			alphaTo[m] ^= mask
		}
		mask <<= 1
	}
	indexOf[alphaTo[m]] = m
	mask >>= 1
	for i := m + 1; i < n; i++ {
		if alphaTo[i - 1] >= mask {
			alphaTo[i] = alphaTo[m] ^ ((alphaTo[i - 1] ^ mask) << 1)
		} else {
			alphaTo[i] = alphaTo[i - 1] << 1
		}
		indexOf[alphaTo[i]] = i
	}
	indexOf[0] = -1

	return
}

func genPoly(t, n, length int, alphaTo, indexOf [1048576]int) (g [548576]int, k, d int, err error) {
	var ii, jj, ll, kaux int
	var aux, nocycles, root, noterms, rdncy int
	var test bool
	cycle, size, min, zeros := make([][]int, 1024), make([]int, 1024), make([]int, 1024), make([]int, 1024)
	for i := 0; i < 1024; i++ {
		cycle[i] = make([]int, 21)
	}

	cycle[0][0] = 0
	size[0] = 1
	cycle[1][0] = 1
	size[1] = 1
	jj = 1
	/*if m > 9 {
		fmt.Printf("Computing cycle sets modulo %d\n", n)
	}*/
	for do := false; !do || ll < n - 1; do = true {
		// Generate the jj-th cycle set
		ii = 0
		for do := false; !do || aux != cycle[jj][0]; do = true {
			ii++
			cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % n
			size[jj]++
			aux = (cycle[jj][ii] * 2) % n
		}
		// Next cycle set representative
		ll = 0
		for do := false; !do || (test && ll < n - 1); do = true {
			ll++
			test = false
			for ii := 1; ii <= jj && !test; ii++ {
				// Examine previous cycle sets
				for kaux = 0; kaux < size[ii] && !test; kaux++ {
					if ll == cycle[ii][kaux] {
						test = true
					}
				}
			}
		}
		if !test {
			jj++ // Next cycle set index
			cycle[jj][0] = ll
			size[jj] = 1
		}
	}
	nocycles = jj // The number of cycle sets modulo n

	/*fmt.Printf("Enter the correcting capability, t: ")
	fmt.Scanf("%d", &t)*/

	d = 2 * t + 1

	// Search for roots 1, 2, ..., d - 1 in cycle sets
	kaux = 0
	rdncy = 0
	for ii = 1; ii <= nocycles; ii++ {
		min[kaux] = 0
		test = false
		for jj = 0; jj < size[ii] && !test; jj++ {
			for root = 1; root < d && !test; root++ {
				if root == cycle[ii][jj] {
					test = true
					min[kaux] = ii
				}
			}
		}
		if min[kaux] != 0 {
			rdncy += size[min[kaux]]
			kaux++
		}
	}
	noterms = kaux
	kaux = 1
	for ii = 0; ii < noterms; ii++ {
		for jj = 0; jj < size[min[ii]]; jj++ {
			zeros[kaux] = cycle[min[ii]][jj]
			kaux++
		}
	}

	k = length - rdncy

	/*if k < 0 {
		fmt.Println("The parameters are invalid!")
	}*/

	fmt.Printf("This is a (%d, %d, %d) binary BCH code.\n", length, k, d)

	// Compute the generator polynomial
	g[0] = alphaTo[zeros[1]]
	g[1] = 1 // g(x) = (x + zeros[1]) initially
	for ii = 2; ii <= rdncy; ii++ {
		g[ii] = 1
		for jj = ii - 1; jj > 0; jj-- {
			if g[jj] != 0 {
				g[jj] = g[jj - 1] ^ alphaTo[(indexOf[g[jj]] + zeros[ii]) % n]
			} else {
				g[jj] = g[jj - 1]
			}
		}
		g[0] = alphaTo[(indexOf[g[0]] + zeros[ii]) % n]
	}
	fmt.Printf("Generator polynomial:\ng(x) = ")
	for ii = 0; ii <= rdncy; ii++ {
		fmt.Printf("%d", g[ii])
		if ii != 0 && ii % 50 == 0 {
			fmt.Printf("\n")
		}
	}
	fmt.Printf("\n")

	return
}

// Compute redundancy bb[], the coefficients of b(x). The redundancy
// polynomial b(x) is the remainder after dividing x^(length-k)*data(x)
// by the generator polynomial g(x).
func encodeBCH(length, k int, g [548576]int, data []int) (bb [548576]int, err error) {
	var i, j, feedback int

	for i = 0; i < length - k; i++ {
		bb[i] = 0
	}
	for i = k - 1; i >= 0; i-- {
		feedback = data[i] ^ bb[length - k - 1]
		if feedback != 0 {
			for j = length - k - 1; j >0; j-- {
				if g[j] != 0 {
					bb[j] = bb[j - 1] ^ feedback
				} else {
					bb[j] = bb[j - 1]
				}
			}
			bb[0] = g[0] & feedback // TODO: was &&
		} else {
			for j = length - k - 1; j > 0; j-- {
				bb[j] = bb[j - 1]
			}
			bb[0] = 0
		}
	}

	return
}

func decodeBCH(length, t, n int, alphaTo, indexOf [1048576]int, recd []int) ([]int, error) {
	var i, j, u, q, t2 int
	var count, synError = 0, 0
	elp, d, l, uLu, s := make([][]int, 1026), make([]int, 1026), make([]int, 1026), make([]int, 1026), make([]int, 1025)
	for i := 0; i < 1024; i++ {
		elp[i] = make([]int, 1024)
	}
	root, loc, /*err,*/ reg := make([]int, 200), make([]int, 200)/*, make([]int, 1024)*/, make([]int, 201)

	t2 = 2 * t

	// First, form the syndromes
	//fmt.Printf("S(x) = ")
	for i = 1; i <= t2; i++ {
		s[i] = 0
		for j = 0; j < length; j++ {
			if recd[j] != 0 {
				s[i] ^= alphaTo[(i * j) % n]
			}
		}
		if s[i] != 0 {
			synError = 1 // Set the error flag if the syndrome is non-zero
			// TODO: If just doing error detection, can stop here
		}

		// Convert the syndrome from polynomial to index form
		s[i] = indexOf[s[i]]
		//fmt.Printf("%3d ", s[i])
	}
	//fmt.Printf("\n")

	// If there are errors, try to correct them
	if synError != 0 {
		// Compute the error location polynomial via the Berlekamp
		// iterative algorithm. Following the terminology of Lin and
		// Costello's book :   d[u] is the 'mu'th discrepancy, where
		// u='mu'+1 and 'mu' (the Greek letter!) is the step number
		// ranging from -1 to 2*t (see L&C),  l[u] is the degree of
		// the elp at that step, and u_l[u] is the difference between
		// the step number and the degree of the elp.

		// Initialize table entries
		d[0] = 0 // Index form
		d[1] = s[1] // Index form
		elp[0][0] = 0 // Index form
		elp[1][0] = 1 // Polynomial form
		for i = 1; i < t2; i++ {
			elp[0][i] = -1 // Index form
			elp[1][i] = 0 // Polynomial form
		}
		l[0] = 0
		l[1] = 0
		uLu[0] = -1
		uLu[1] = 0
		u = 0

		for do := false; !do || (u < t2 && l[u + 1] <= t); do = true {
			u++
			if d[u] == -1 {
				l[u + 1] = l[u]
				for i = 0; i <= l[u]; i++ {
					elp[u + 1][i] = elp[u][i]
					elp[u][i] = indexOf[elp[u][i]]
				}
			} else {
				// Search for words with greatest u_lu[q] for which d[q] != 0
				q = u - 1
				for d[q] == -1 && q > 0 {
					q--
				}
				// Have found first non-zero d[q]
				if q > 0 {
					j = q
					for do := false; !do || j > 0; do = true {
						j--
						if d[j] != -1 && uLu[q] < uLu[j] {
							q = j
						}
					}
				}

				// Have found q such that d[u] != 0 and u_lu[q] is maximum
				// Store the degree of new elp polynomial
				if l[u] > l[q] + u - q {
					l[u + 1] = l[u]
				} else {
					l[u + 1] = l[q] + u - q
				}

				// Form new elp(x)
				for i = 0; i < t2; i++ {
					elp[u + 1][i] = 0
				}
				for i = 0; i <= l[q]; i++ {
					if elp[q][i] != -1 {
						elp[u + 1][i + u - q] = alphaTo[(d[u] + n - d[q] + elp[q][i]) % n]
					}
				}
				for i = 0; i <= l[u]; i++ {
					elp[u + 1][i] ^= elp[u][i]
					elp[u][i] = indexOf[elp[u][i]]
				}
			}
			uLu[u + 1] = u - l[u + 1]

			// Form the (u + 1)th discrepancy
			if u < t2 {
				// No discrepancy computed on the last iteration
				if s[u + 1] != -1 {
					d[u + 1] = alphaTo[s[u + 1]]
				} else {
					d[u + 1] = 0
				}
				for i = 1; i <= l[u + 1]; i++ {
					if s[u + 1 - i] != -1 && elp[u + 1][i] != 0 {
						d[u + 1] ^= alphaTo[(s[u + 1 - i] + indexOf[elp[u + 1][i]]) % n]
					}
				}
				// Put d[u + 1] into index form
				d[u + 1] = indexOf[d[u + 1]]
			}
		}

		u++
		if l[u] <= t { // Can correct errors
			// Put elp into index form
			for i = 0; i <= l[u]; i++ {
				elp[u][i] = indexOf[elp[u][i]]
			}

			/*(fmt.Printf("sigma(x) = ")
			for i = 0; i <= l[u]; i++ {
				fmt.Printf("%3d ", elp[u][i])
			}
			fmt.Printf("\n")
			fmt.Printf("Roots: ")*/

			// Chien search: find the roots of the error location polynomial
			for i = 1; i <= l[u]; i++ {
				reg[i] = elp[u][i]
			}
			count = 0
			for i = 1; i <= n; i++ {
				q = 1
				for j = 1; j <= l[u]; j++ {
					if reg[j] != -1 {
						reg[j] = (reg[j] + j) % n
						q ^= alphaTo[reg[j]]
					}
				}
				// Store the root and error location number indices
				if q == 0 {
					root[count] = i
					loc[count] = n - i
					count++
					//fmt.Printf("%3d ", n - i)
				}
			}
			//fmt.Printf("\n")
			if count == l[u] {
				// Number of roots = degree of elp, hence <= t errors
				for i = 0; i < l[u]; i++ {
					recd[loc[i]] ^= 1
				}
			} else {
				fmt.Println("Incomplete decoding - errors detected.")
			}
		}
	}

	return recd[:], nil
}

func Encode(codeLength, correctableErrors int, buf *[]int) ([]int, EncodingConfig) {
	m := int(math.Log2(float64(codeLength))) + 1

	n, p, _ := readP(m)
	alphaTo, indexOf, _ := generateGF(m, n, p)
	g, k, _, _ := genPoly(correctableErrors, n, codeLength, alphaTo, indexOf)

	data := make([]int, k)
	for i := 0; i < k; i++ {
		data[i] = (*buf)[i]
	}

	bb, _ := encodeBCH(codeLength, k, g, data)

	// Load the data into the code value (ecc bits followed by original data)
	recd := make([]int, codeLength)
	for i := 0; i < codeLength - k; i++ {
		recd[i] = bb[i]
	}
	for i := 0; i < k; i++ {
		recd[i + codeLength - k] = data[i]
	}

	for i := 0; i < codeLength; i++ {
		fmt.Printf("%1d", recd[i])
		if i != 0 && i % 50 == 0 {
			fmt.Printf("\n")
		}
	}
	fmt.Printf("\n")

	return recd, EncodingConfig{correctableErrors, n, GField{alphaTo, indexOf}}
}

func Decode(buf *[]int, config EncodingConfig) []int {
	recd, _ := decodeBCH(len(*buf), config.MaxCorrectableErrors, config.N, config.Field.AlphaTo, config.Field.IndexOf, *buf)

	return recd
}