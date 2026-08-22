use std::error::Error;

#[test]
fn test_valid_read() -> Result<(), Box<dyn Error>> {
    use std::path::PathBuf;

    use bigtools::BigWigRead;

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut valid_bigwig = dir.clone();
    valid_bigwig.push("valid.bigWig");

    let mut bwread = BigWigRead::open_file(valid_bigwig).unwrap();

    // Test that chrom tree parsing works
    let chroms = bwread.chroms();
    assert_eq!(chroms.len(), 1);
    // chr17
    assert_eq!(chroms[0].length, 83257441);

    // Test that header/summary is correct
    let summary = bwread.get_summary()?;
    assert_eq!(summary.bases_covered, 137894);
    assert_eq!(summary.max_val, 14254.0);

    // This tests simply reading an interval. Importantly, the provided interval actually splits an interval in two, so this also tests correct splitting on read.
    let first_interval = bwread
        .get_interval("chr17", 0, 59899)?
        .next()
        .unwrap()
        .unwrap();
    assert_eq!(first_interval.start, 59898);
    assert_eq!(first_interval.end, 59899);
    assert_eq!(first_interval.value, 0.06792);

    Ok(())
}

#[test]
fn test_values() -> Result<(), Box<dyn Error>> {
    use std::path::PathBuf;

    use bigtools::BigWigRead;

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut valid_bigwig = dir.clone();
    valid_bigwig.push("valid.bigWig");

    let mut bwread = BigWigRead::open_file(valid_bigwig).unwrap();

    let vals = bwread.values("chr17", 0, 59899)?;
    assert_eq!(vals.len(), 59899);
    assert_eq!(vals[59898], 0.06792);
    Ok(())
}

#[test]
fn test_reduction_values() -> Result<(), Box<dyn Error>> {
    use std::path::PathBuf;

    use bigtools::BigWigRead;

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    let mut valid_bigwig = dir.clone();
    valid_bigwig.push("valid.bigWig");

    let mut bwread = BigWigRead::open_file(valid_bigwig).unwrap();
    let interval = bwread.get_zoom_interval("chr17", 0, 36996442, 10240);
    let x: Vec<_> = interval.unwrap().collect();

    assert_eq!(x.len(), 16);
    Ok(())
}

#[test]
fn test_primary_data_block_layout() -> Result<(), Box<dyn Error>> {
    use std::path::PathBuf;

    use bigtools::{BigBedRead, BigWigRead};

    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push("resources/test");

    for (file_name, is_bigwig) in [("valid.bigWig", true), ("bigGenePred.bb", false)] {
        let path = dir.join(file_name);
        let (blocks, chrom_ids) = if is_bigwig {
            let mut reader = BigWigRead::open_file(path)?;
            let chrom_ids: Vec<_> = reader.chroms().iter().map(|chrom| chrom.id()).collect();
            (reader.data_blocks()?, chrom_ids)
        } else {
            let mut reader = BigBedRead::open_file(path)?;
            let chrom_ids: Vec<_> = reader.chroms().iter().map(|chrom| chrom.id()).collect();
            (reader.data_blocks()?, chrom_ids)
        };

        assert!(!blocks.is_empty(), "{file_name} should contain data blocks");
        assert!(
            blocks.iter().all(|block| block.compressed_size > 0),
            "{file_name} reported an empty compressed block"
        );
        assert!(blocks.windows(2).all(|pair| {
            (pair[0].start_chrom_id, pair[0].start_base)
                <= (pair[1].start_chrom_id, pair[1].start_base)
        }));
        assert!(blocks.iter().all(|block| {
            chrom_ids.contains(&block.start_chrom_id) && chrom_ids.contains(&block.end_chrom_id)
        }));
    }

    Ok(())
}
