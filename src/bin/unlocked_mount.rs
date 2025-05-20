//
// GRSIM script to simulate the effect of RLE ground accelerations
//

use gmt_dos_actors::actorscript;
use gmt_dos_clients::Source; //{Signals, Signal, }
use gmt_dos_clients_fem::{solvers::ExponentialMatrix, DiscreteModalSolver};
use gmt_dos_clients_io::{
    gmt_fem::{
        inputs::OSS00GroundAcc,
        outputs::{OSS00Ground6D, OSSPayloads6D, Pier6D}, //,OSSHardpointForce
    },
    mount::{MountEncoders, MountTorques}, // MountSetPoint,
};

use gmt_dos_clients_mount::Mount;
use gmt_fem::FEM;

use matio_rs::MatFile;
use nalgebra as na;
use std::{env, path::Path};

/*
export FEM_REPO=$HOME/mnt/20250313_0917_zen_30_M1_202110_FSM_202305_Mount_202305_pier_202411_heavy_SIIDOM_movingGround
export MOUNT_MODEL=MOUNT_FDR_1kHz
cargo run --release --bin unlocked_mount
*/

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    unsafe {
        env::set_var(
            "DATA_REPO",
            Path::new(env!("CARGO_MANIFEST_DIR")).join("data"),
        );
    }
    let sim_sampling_frequency = 1000;
    //let n_step = sssha_length;//10000; //
    let fem = FEM::from_env()?;

    for id in 1..=2 {
        print!("RLE #{:02} dataset: ", id);
        // RLE data path
        let sssha_path = Path::new(&env::var("CARGO_MANIFEST_DIR").unwrap())
            .join("data")
            .join(format!("RLE{:02}_1kHz_dt.mat", id));
        let acc_mat: na::DMatrix<f64> =
            MatFile::load(&sssha_path)?.var(format!("RLE{:02}_1kHz_dt", id))?;
        let (_sssha_length, n_sssha_dim) = acc_mat.shape();
        println!("{:?}", acc_mat.shape());

        // DISCRETE STATE-SPACE MODEL (from FEM)
        let state_space = DiscreteModalSolver::<ExponentialMatrix>::from_fem(fem.clone())
            .sampling(sim_sampling_frequency as f64)
            .proportional_damping(2. / 100.)
            // .max_eigen_frequency(95f64)
            .including_mount()
            .ins::<OSS00GroundAcc>()
            .outs::<Pier6D>()
            .outs::<OSS00Ground6D>()
            .outs::<OSSPayloads6D>()
            //.use_static_gain_compensation()
            .build()?;
        println!("{state_space}");

        // FEM
        let fem = state_space;
        // MOUNT CONTROL
        let mount = Mount::new();

        // RLE-SSSHA DATA
        let sssha_source = Source::new(acc_mat.transpose().as_slice().to_vec(), n_sssha_dim);

        // DSL
        actorscript! {
        #[labels(fem = "Telescope\nStructure", sssha_source = "SSSHA Acc")] //mount = "Mount\nControl", 
        1: sssha_source[OSS00GroundAcc] -> fem
        // Mount control loop
        //1: setpoint[MountSetPoint] -> mount[MountTorques] -> fem[MountEncoders]! -> mount
        //
        //1: fem[MountEncoders]! -> mount[MountTorques] -> fem        
        // LOG
        1: fem[OSSPayloads6D]${162}
        1: sssha_source[OSS00GroundAcc]${3}
        1: fem[Pier6D]${12}
        1: fem[OSS00Ground6D]${6}
        1: fem[MountEncoders]${14}
        }

        let fem_id = Path::new(env!("FEM_REPO"))
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .splitn(3, "_")
            .take(2)
            .collect::<Vec<_>>()
            .join("_");
        model_logging_1
            .lock()
            .await
            .to_parquet(format!("model-{}-RLE{:02}_mntOL.parquet",fem_id, id))?;
    }
    Ok(())
}
